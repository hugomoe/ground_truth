#define ZOOM 10
#include "fft_zoom.h"


int apply_homo_ground_truth(float *img,float *img_f,int w,int h,int w_f,int h_f,double H[3][3]){	
	
		double fzoom=ZOOM;
		double HH[3][3];
		HH[0][0]=H[0][0];
		HH[0][1]=H[0][1];
		HH[0][2]=H[0][2]*fzoom;
		HH[1][0]=H[1][0];
		HH[1][1]=H[1][1];
		HH[1][2]=H[1][2]*fzoom;
		HH[2][0]=H[2][0]/fzoom;
		HH[2][1]=H[2][1]/fzoom;
		HH[2][2]=H[2][2];


	
	float *img_aux = malloc(3*sizeof(float)*w*h*ZOOM*ZOOM);
	float *img_aux2 = malloc(3*sizeof(float)*w_f*h_f*ZOOM*ZOOM);
	for(int i=0;i<w*h*ZOOM*ZOOM*3;i++){img_aux[i]=0;}
	
	
	
	/*for(int l=0;l<3;l++){
	for(int i=0;i<w;i++){
		for(int j=0;j<h;j++){
			for(int u=0;u<ZOOM;u++){
				for(int v=0;v<ZOOM;v++){
					float x=u/fzoom;
					float y=v/fzoom;
					int id,jd;
					if(i==w-1){id=0;}else{id=i+1;}
					if(j==h-1){jd=0;}else{jd=j+1;}
					img_aux[3*(i*ZOOM+u+(j*ZOOM+v)*w*ZOOM)+l]=(1-x)*(1-y)*(img[(i+j*w)*3+l])+(1-x)*y*(img[(i+jd*w)*3+l])
					+x*(1-y)*(img[(id+j*w)*3+l])+x*y*(img[(id+jd*w)*3+l]);
				}
			}
		}
	}
	}*/
	
	int w_zoom = ZOOM*w; int h_zoom = ZOOM*h;
	
	zoom(img,w,h,3,w_zoom,h_zoom,img_aux);
	
	
	
	
	interpolator_t EVAL = bilinear_interpolation_at; 
	extrapolator_t OUT = getsample_cons;
	
	
	for(int l=0;l<3;l++)
	for (int j = 0; j < ZOOM*h_f; j++)
	for (int i = 0; i < ZOOM*w_f; i++)
	{
		double p[2] ={i,j};
		
		apply_homography(p, HH, p);
		p[0] = (p[0] - 0.5) * ZOOM * w / (ZOOM * w - 1.0);
		p[1] = (p[1] - 0.5) * ZOOM * h / (ZOOM * h - 1.0);
			int idx = 3*(ZOOM * w_f * j + i)+l;
			img_aux2[idx] = EVAL(img_aux, ZOOM*w, ZOOM*h, 3, p[0], p[1], l, OUT);
	}
	
	int taps = 2*ZOOM;
	double sigma = 0.8 * ZOOM;
	double *gauss = malloc(pow((2*taps+1),2)*sizeof(float));
	double tot = 0;
	for(int i=-taps;i<=taps;i++){
		for(int j=-taps;i<=taps;i++){
			tot += (gauss[i+taps+(j+taps)*(2*taps+1)] = exp(-(pow(i,2)+pow(j,2))/(2*pow(sigma,2))));
		}
	}	
	for(int u=0;u<pow(2*taps+1,2);u++){gauss[u]=gauss[u]/tot;}

	
	for(int l=0;l<3;l++){
		for (int j = 0; j < h_f; j++){
			for (int i = 0; i < w_f; i++){
				float v=0;
				for(int i2=-taps;i2<=taps;i2++){
					for(int j2=-taps;j2<=taps;j2++){
						int i1 = good_modulus(i*ZOOM+i2,w_f*ZOOM);
						int j1 = good_modulus(j*ZOOM+j2,h_f*ZOOM);
						v += gauss[i2+taps+(j2+taps)*(2*taps+1)]*img_aux2[(i1+j1*ZOOM*w_f)*3+l];
						int idx = l + 3 * (w_f * j + i);
						img_f[idx] = v;
					}
				}
			}
		}
	}
	
	free(gauss);
	free(img_aux);
	free(img_aux2);
	
	return 0;
}
