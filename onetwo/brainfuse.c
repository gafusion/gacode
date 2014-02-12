#include <stdio.h>
#include "doublefann.h"

int run_net_on_data_(int * num_data, int * num_input, int * num_output, double * input, double * output, char *path, int *debug){
  struct fann *ann;
  fann_type *calc_in, *calc_out;
  int i, j, k;

  calc_in  = (fann_type *)malloc(*num_input  * sizeof(fann_type));
  printf("n=%d in=%d out=%d prec=%d\n",*num_data,*num_input,*num_output,sizeof(fann_type));

  //////////////////////////////
  if (debug) printf("- network from %s\n",path);
  //////////////////////////////
  ann = fann_create_from_file(path);
  if(!ann){
      printf("Error creating ann --- ABORTING.\n");
      return -1;
  }

  //////////////////////////////
  if (debug) printf("- calculate\n");
  //////////////////////////////
  for(i = 0; i < *num_data; i++){
    //load
    for(j = 0; j < *num_input; j++){
      calc_in[j]=(fann_type)input[i + j * *num_data];
    }
    //scale in
    fann_scale_input(ann,calc_in);
    //compute
    calc_out = fann_run(ann, calc_in);
    //scale out
    fann_descale_output(ann,calc_out);
    //save
    for(j = 0; j < *num_output; j++){
      output[i + j * *num_data]=(double)calc_out[j];
    }
  }

  //////////////////////////////
  if(debug>1){
    for(i = 0; i < *num_data; i++){
      for( j = 0; j < *num_output; j++){
	printf("%f ",output[i + j * *num_data]);
      }
      printf("\n");
    }
  }
  //////////////////////////////

  //////////////////////////////
  if (debug) printf("- cleanup\n");
  //////////////////////////////
  fann_destroy(ann);
  free(calc_in);

  return 0;
}

