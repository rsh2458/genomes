/****************************************************************
 ****
 ****               Final Project - OpenMP
 ****
 ****  Comparing a Set of Genomes to a Reference Genome
 ****
 ****                     Rowan Hart
 ****
 ****************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <dirent.h>


// *****************************************
//  CONFIGURE - - Number of Query Genomes
// *****************************************
#define number_of_query_genomes 115
// *****************************************


// Find the size of the reference genome
long int findSize(char file_name[])
{
    FILE* fp = fopen(file_name, "r");
    if (fp == NULL) {
        printf("File Not Found!\n");
        return -1;
    }
    fseek(fp, 0L, SEEK_END);
    long int res = ftell(fp);
    fclose(fp);
    return res;
}

// Creates an array of paths from the input file
void find_paths(char* fileName, char pathways[number_of_query_genomes][256]){
    FILE* file = fopen(fileName, "r"); /* should check the result */
    char line[256];
    int init = 0;
    while (fgets(line, sizeof(line), file)) {
        line[strcspn(line, "\n")] = 0;
        strcpy(pathways[init], line);
        init ++;
    }
    fclose(file);
}

// Create a char array for the reference genome + Calculating the G,C,T,A percentages for refernce genome
void nucleotide_array(char array[], char* filename, long int* Gcount, long int* Ccount, long int* Acount, long int* Tcount){
    printf("Creating Reference Array... \n");
    FILE *fp;
    char c;
    Ccount=0;
    Gcount=0;
    Acount=0;
    Tcount=0;
    long int init = 0;
    fp = fopen(filename, "r");
    if (fp == NULL){
        printf("Could not open file %s\n",filename);
        return;
    }
    while(1){
        c = fgetc(fp); // read the file
        if (c == EOF) break;
        // G count
        if (c == 'G') Gcount++;

        // A count
        if (c == 'A') Acount++;

        // C count
        if (c == 'C') Ccount++;

        // T count
        if (c == 'T') Tcount++;
        array[init] = c;
        init += 1;
    }
    fclose(fp);
    return;
}


int main() {

// Initialize the Parallelism 

    printf("\n\nHello from master thread\n\n");
    int num_threads;

    #pragma omp parallel
    {
    printf("This is thread %d\n", omp_get_thread_num());
    #pragma omp single
    {
    num_threads = omp_get_num_threads();
    }
    }

// Display that the program is running
    printf("\n");
    printf("Rowan Hart's Final Project is running...\n");
   // printf("\n");

// *****************************************
//   CONFIGURE -- Pathways to Genomes
// *****************************************

// Pathway to reference genome
    char* ref_filename = "/home1/07947/rhart/final_project/reference_genomes/Salmonella";

// Pathway to list of query genomes and number of genomes
    char* list_of_query_genomes = "/home1/07947/rhart/final_project/list_of_query_genomes.list";

// *****************************************
// *****************************************

// Count the number of nucleotides in the reference genome
    long int ref_length = findSize(ref_filename);

// Create a char array for the reference genome
    long int ref_G, ref_A, ref_C, ref_T;
    char* ref_genome;
    ref_genome = (char*)malloc(ref_length * sizeof(char)); // Malloc array
    ref_genome[ref_length] = EOF;
    nucleotide_array(ref_genome, ref_filename, &ref_G, &ref_C, &ref_A, &ref_T);

// Read the List of Query Genome Paths
    char query_paths[number_of_query_genomes][256];
    find_paths(list_of_query_genomes, query_paths);

    printf("Comparing Genomes...\n");

    double f1 = omp_get_wtime();
// Initialize Characteristics of Each Genome
    float ANI[number_of_query_genomes];
    float G_count[number_of_query_genomes];
    float C_count[number_of_query_genomes];
    float A_count[number_of_query_genomes];
    float T_count[number_of_query_genomes];
    long int nt_size[number_of_query_genomes];

// Initialize Variables
    long int Gcount, Ccount, Acount, Tcount, match, j;
    char* curr_query;
    char nuc;
    FILE *fp;

// Initialize Genome counter
    int genome_count[num_threads];
    #pragma omp parallel
    {
        genome_count[omp_get_thread_num()] = 0;
    }

// Compare Genomes
    
    #pragma omp parallel for private(j, curr_query, nuc, fp, Acount, Gcount, Tcount, Ccount, match) schedule(dynamic, 1)
    for (int i = 0; i<number_of_query_genomes; i++){
        // Initialize Variables to 0
        j=0;
        match=0;
        Gcount=0;
        Acount=0;
        Ccount=0;
        Tcount=0;

        // Add to thread genome total
        genome_count[omp_get_thread_num()] += 1;

        // Declare what genome a thread is working on
        curr_query = query_paths[i];

        // Calculate Length of Query Genome
        nt_size[i] = findSize(curr_query);

        // Compare Query Genome to Reference Genome

        fp = fopen(curr_query, "r");
        if (fp == NULL){
            printf("Could not open file %s\n",curr_query);
        }
        else{
            while(1){
                nuc = fgetc(fp); // read a character from the file
                if (nuc == EOF || ref_genome[j] == EOF) break; // Check to see if the file has ended

                // Compare Nucleotides
                if (nuc == ref_genome[j]) match++;
                // G count
                if (nuc == 'G') Gcount++;
                // A count
                if (nuc == 'A') Acount++;
                // C count
                if (nuc == 'C') Ccount++;
                // T count
                if (nuc == 'T') Tcount++;
                j++;
            }
        }
        fclose(fp);
        
        // Return Values
        G_count[i] = (float)Gcount / (float)j;
        C_count[i] = (float)Ccount / (float)j;
        A_count[i] = (float)Acount / (float)j;
        T_count[i] = (float)Tcount / (float)j;
        ANI[i] = (float)match / (float)j;
    }
    double f2 = omp_get_wtime();


// Print Results


    printf("\n");
    printf("----------------------------------------------------\n");
    printf("Comparing Genomes - Parallel Computing - Rowan Hart \n");
    printf("----------------------------------------------------\n");

    printf("\n");
    printf("%-35s\n","OpenMP");
    printf("--------------------------------------\n");
    printf("%-22s:: %10d\n","Num Threads", num_threads);
    printf("--------------------------------------\n");
    printf("\n");

    printf("\n");
    printf("------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    printf("%-110s:: %-12s :: %-20s :: %-12s :: %-12s :: %-12s :: %12s ::\n", "Reference Genome", "ANI (%)", "Nucleotides (bytes)", "G Count (%)", "C Count (%)", "A Count (%)", "T Count (%)" );
    printf("------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n\n");
    printf("%-110s:: %-12.2f :: %-20d :: %-12.2f :: %-12.2f :: %-12.2f :: %12.2f ::\n\n", ref_filename, 100, ref_length, ref_G * 100, ref_C * 100, ref_A * 100, ref_T * 100);
    printf("------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    printf("%-110s:: %-12s :: %-20s :: %-12s :: %-12s :: %-12s :: %12s ::\n", "Query Genomes (filename pathway)", "ANI (%)", "Nucleotides (bytes)", "G Count (%)", "C Count (%)", "A Count (%)", "T Count (%)" );
    printf("------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n\n");
    for (int i = 0; i<number_of_query_genomes; i++){
    printf("%-110s:: %-12.2f :: %-20d :: %-12.2f :: %-12.2f :: %-12.2f :: %12.2f ::\n", query_paths[i], ANI[i] * 100, nt_size[i], G_count[i] * 100, C_count[i] * 100, A_count[i] * 100, T_count[i] * 100);
    }
    printf("------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n\n");

    printf("\n");
    printf("%-22s:: %10s\n","Parallel Action", "time/s");
    printf("-----------------------------------\n");
    printf("%-22s:: %9.3f\n", "Comparing Genomes", (f2-f1));
    printf("-----------------------------------\n");
    printf("\n");

    printf("\n");
    printf("%-22s:: %10s\n","Thread Number", "Genomes Compared");
    printf("-----------------------------------\n");
    #pragma omp parallel
    {
    printf("%-22d:: %10d\n", omp_get_thread_num(), genome_count[omp_get_thread_num()]);
    }
    printf("-----------------------------------\n");
    printf("\n");

  return 0;
}
