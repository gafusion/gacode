void func(float ri, float phi, float theta, unsigned short *s, float *x)
{
  complex<float> f_c;
  complex<float> cs;
  complex<float> IC(0.0,1.0);
  float sum;

  int i_loc;

  //---------------------------------------------------
  //x[0] = (phi+q[i_r]*theta)/ly-0.5;
  //x[1] = (ri-r[0])/(r[n_r-1]-r[0])-0.5;
  //x[2] = theta/ly;
  //---------------------------------------------------

  sum = 0.0;

  for (int i_n=i_n_0; i_n<n_n; i_n++) {

    i_loc = i_r*(n_theta_plot)+ 
      i_field*(n_r*n_theta_plot) +
      i_n*(n_field*n_r*n_theta_plot) + 
      i_time*(n_n*n_field*n_r*n_theta_plot);

    f_c = (f_real[i_loc]+IC*f_imag[i_loc]); 

    cs = exp(-(float)(n[i_n])*IC*(phi-q[i_r]*theta));

   if (i_n > 0) {
      sum = sum+real(cs*f_c);
    } else {
      sum = sum+0.5*real(cs*f_c);
    }

  }

  *s = 128+(unsigned short) (function_magnification*sum);
}
