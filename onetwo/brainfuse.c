#include <stdio.h>
#include "fann.h"

int run_net_on_data_(int * num_data, int * num_input, int * num_output, fann_type * input, fann_type * output){
  struct fann *ann;
  fann_type *calc_in, *calc_out;
  int i, j, k;

  calc_in  = (fann_type *)malloc(*num_input  * sizeof(fann_type));
  printf("%d %d %d\n",*num_data,*num_input,*num_output);

  //////////////////////////////
  printf("- network\n");
  //////////////////////////////
  ann = fann_create_from_file("brainfuse.net");
  if(!ann){
      printf("Error creating ann --- ABORTING.\n");
      return -1;
  }

  //////////////////////////////
  printf("- calculate\n");
  //////////////////////////////
  for(i = 0; i < *num_data; i++){
    for(j = 0; j < *num_input; j++){
      calc_in[j]=input[i + j * *num_data];
    }
    calc_out = fann_run(ann, calc_in);
    for(j = 0; j < *num_output; j++){
      output[i + j * *num_data]=calc_out[j];
    }
  }

  //for(i = 0; i < *num_data; i++){
  //  printf("  (%f, %f) -> %f\n", input[i], input[i+*num_data], output[i]);
  //}

  //////////////////////////////
  printf("- cleanup\n");
  //////////////////////////////
  fann_destroy(ann);
  free(calc_in);

  return 0;
}

/*

int hello_(char *cc, int ll){
  cc[ll--] = '\0';  // NULL terminate the string
  printf("hello %s\n",cc);
  return 0;
}

int main1(){
  int ret;
  unsigned int i;
  int num_data, num_input, num_output;

  //////////////////////////////
  printf("DATA\n");
  //////////////////////////////
  num_data=4;
  num_input=2;
  num_output=1;
  fann_type ** input  = (fann_type **) calloc(num_data, sizeof(fann_type *));
  fann_type ** output = (fann_type **) calloc(num_data, sizeof(fann_type *));
  for (i = 0; i < num_data; i++){
    input[i]  = (fann_type *)malloc(num_input  * sizeof(fann_type));
    output[i] = (fann_type *)malloc(num_output * sizeof(fann_type));
  }

  input[0][0]=-1;
  input[0][1]=-1;
  output[0][0]=-1;

  input[1][0]=-1;
  input[1][1]=1;
  output[1][0]=1;

  input[2][0]=1;
  input[2][1]=-1;
  output[2][0]=1;

  input[3][0]=1;
  input[3][1]=1;
  output[3][0]=-1;


  //////////////////////////////
  printf("CALL\n");
  //////////////////////////////
  run_net_on_data_(num_data,num_input,num_output,input,output);

  //////////////////////////////
  printf("CHECK\n");
  //////////////////////////////
  for(i = 0; i < num_data; i++){
    printf("XOR test (%f, %f) -> %f\n", input[i][0], input[i][1], output[i][0]);
  }

  return ret;
}

*/
