/*******************************************************************************
% Code originally by Philippe Weinzaepfel, 2015, 
% for motion boundaries detection
% Licensed under the MSR-LA Full Rights License [see license.txt]
*******************************************************************************/

#include <mex.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include <stdint.h>

#define MIN(a,b)  (((a)<(b)) ? (a) : (b))
#define MAX(a,b)  (((a)>(b)) ? (a) : (b))
#define NEWA(type,n) (type*)malloc(sizeof(type)*(n))

/* does not seems to be included on windows */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static inline float pow2( float f ) {
  return f*f;
}

template<typename type_in_t, typename type_out_t>
void _smooth_3_horiz(int tx, int ty, const int w_center, const int w_side, type_in_t* pixels, type_out_t* _res) {
  int y;
  int sum_w = 2*w_side + w_center;
  if(!sum_w)  sum_w=1;
  #if defined(USE_OPENMP)
  #pragma omp parallel for
  #endif
  for(y=0; y<ty; y++) {
    int x,pos = y*tx;
    type_out_t* res = _res + pos;
    *res++ = ( (w_center+w_side)*pixels[0+pos] + w_side*pixels[1+pos])/sum_w;
    for(x=1; x<tx-1; x++)
      *res++ = (w_side*pixels[x+1+pos] + w_center*pixels[x+pos] + w_side*pixels[x-1+pos])/sum_w;
    *res++ = ( (w_center+w_side)*pixels[x+pos] + w_side*pixels[x-1+pos])/sum_w;
  }
}


template<typename type_in_t, typename type_out_t>
void _smooth_5_horiz(int tx, int ty, 
                     const int w_center, const int w_side1, const int w_side2, 
                     type_in_t* pixels, type_out_t* _res) {
  int y;
  int sum_w = 2*(w_side1 + w_side2) + w_center;
  if(!sum_w)  sum_w=1;
  #if defined(USE_OPENMP)
  #pragma omp parallel for
  #endif
  for(y=0; y<ty; y++) {
    int x,pos = y*tx;
    type_out_t* res = _res + pos;
      x=0;
      *res++ = ( 
                w_side2 * pixels[x  +pos] +
                w_side1 * pixels[x  +pos] +
                w_center* pixels[x  +pos] +
                w_side1 * pixels[x+1+pos] +
                w_side2 * pixels[x+2+pos] ) / sum_w;
      x++;
      *res++ = ( 
                w_side2 * pixels[x-1+pos] +
                w_side1 * pixels[x-1+pos] +
                w_center* pixels[x  +pos] +
                w_side1 * pixels[x+1+pos] +
                w_side2 * pixels[x+2+pos] ) / sum_w;
    
    for(x=2; x<tx-2; x++)
      *res++ = ( 
                w_side2 * pixels[x-2+pos] +
                w_side1 * pixels[x-1+pos] +
                w_center* pixels[x  +pos] +
                w_side1 * pixels[x+1+pos] +
                w_side2 * pixels[x+2+pos] ) / sum_w;
    
      *res++ = ( 
                w_side2 * pixels[x-2+pos] +
                w_side1 * pixels[x-1+pos] +
                w_center* pixels[x  +pos] +
                w_side1 * pixels[x+1+pos] +
                w_side2 * pixels[x+1+pos] ) / sum_w;
      x++;
      *res++ = ( 
                w_side2 * pixels[x-2+pos] +
                w_side1 * pixels[x-1+pos] +
                w_center* pixels[x  +pos] +
                w_side1 * pixels[x  +pos] +
                w_side2 * pixels[x  +pos] ) / sum_w;
  }
}

template<typename type_in_t, typename type_out_t>
void _smooth_7_horiz(int tx, int ty, const int w_center, const int w_side1, const int w_side2, const int w_side3, 
                    type_in_t* pixels, type_out_t* _res) {
  int y;
  int sum_w = 2*(w_side1 + w_side2 + w_side3) + w_center;
  if(!sum_w)  sum_w=1;
  #if defined(USE_OPENMP)
  #pragma omp parallel for
  #endif
  for(y=0; y<ty; y++) {
    int x,pos = y*tx;
    type_out_t* res = _res + pos;
      x=0;
      *res++ = ( 
                w_side3 * pixels[x  +pos] +
                w_side2 * pixels[x  +pos] +
                w_side1 * pixels[x  +pos] +
                w_center* pixels[x  +pos] +
                w_side1 * pixels[x+1+pos] +
                w_side2 * pixels[x+2+pos] +
                w_side3 * pixels[x+3+pos] ) / sum_w;
      x++;
      *res++ = ( 
                w_side3 * pixels[x-1+pos] +
                w_side2 * pixels[x-1+pos] +
                w_side1 * pixels[x-1+pos] +
                w_center* pixels[x  +pos] +
                w_side1 * pixels[x+1+pos] +
                w_side2 * pixels[x+2+pos] +
                w_side3 * pixels[x+3+pos] ) / sum_w;
      x++;
      *res++ = ( 
                w_side3 * pixels[x-2+pos] +
                w_side2 * pixels[x-2+pos] +
                w_side1 * pixels[x-1+pos] +
                w_center* pixels[x  +pos] +
                w_side1 * pixels[x+1+pos] +
                w_side2 * pixels[x+2+pos] +
                w_side3 * pixels[x+3+pos] ) / sum_w;
    
    for(x=3; x<tx-3; x++)
      *res++ = ( 
                w_side3 * pixels[x-3+pos] +
                w_side2 * pixels[x-2+pos] +
                w_side1 * pixels[x-1+pos] +
                w_center* pixels[x  +pos] +
                w_side1 * pixels[x+1+pos] +
                w_side2 * pixels[x+2+pos] +
                w_side3 * pixels[x+3+pos] ) / sum_w;
    
      *res++ = ( 
                w_side3 * pixels[x-3+pos] +
                w_side2 * pixels[x-2+pos] +
                w_side1 * pixels[x-1+pos] +
                w_center* pixels[x  +pos] +
                w_side1 * pixels[x+1+pos] +
                w_side2 * pixels[x+2+pos] +
                w_side3 * pixels[x+2+pos] ) / sum_w;
      x++;
      *res++ = ( 
                w_side3 * pixels[x-3+pos] +
                w_side2 * pixels[x-2+pos] +
                w_side1 * pixels[x-1+pos] +
                w_center* pixels[x  +pos] +
                w_side1 * pixels[x+1+pos] +
                w_side2 * pixels[x+1+pos] +
                w_side3 * pixels[x+1+pos] ) / sum_w;
      x++;
      *res++ = ( 
                w_side3 * pixels[x-3+pos] +
                w_side2 * pixels[x-2+pos] +
                w_side1 * pixels[x-1+pos] +
                w_center* pixels[x  +pos] +
                w_side1 * pixels[x  +pos] +
                w_side2 * pixels[x  +pos] +
                w_side3 * pixels[x  +pos] ) / sum_w;
  }
}

template<typename type_in_t, typename type_out_t>
void _smooth_3_vert(int tx, int ty, const int w_center, const int w_side, type_in_t* pixels, type_out_t* res) {
  int x,y,pos=0;
  int sum_w = 2*w_side + w_center;
  if(!sum_w)  sum_w=1;
  for(x=0; x<tx; x++,pos++)
    res[pos] = (  (w_center+w_side)*pixels[pos] + w_side*pixels[pos+tx])/sum_w;
  #if defined(USE_OPENMP)
  #pragma omp parallel for
  #endif
  for(y=1; y<ty-1; y++) {
    int x,pos = y*tx;
    for(x=0; x<tx; x++,pos++)
      res[pos] = ( w_side*pixels[pos+tx] + w_center*pixels[pos] + w_side*pixels[pos-tx])/sum_w;
  }
  pos = (ty-1)*tx;
  for(x=0; x<tx; x++,pos++)
    res[pos] = ( (w_center+w_side)*pixels[pos] + w_side*pixels[pos-tx])/sum_w;
}

template<typename type_in_t, typename type_out_t>
void _smooth_5_vert(int tx, int ty, const int w_center, const int w_side1, const int w_side2, type_in_t* pixels, type_out_t* res) {
  int x,y,pos=0;
  int sum_w = 2*(w_side1 + w_side2) + w_center;
  if(!sum_w)  sum_w=1;
  const int tx1=tx,tx2=2*tx;
  for(x=0; x<tx; x++,pos++)
      res[pos] = ( 
                  w_side2 * pixels[pos] + 
                  w_side1 * pixels[pos] + 
                  w_center* pixels[pos] + 
                  w_side1 * pixels[pos+tx1] + 
                  w_side2 * pixels[pos+tx2] 
                  )/sum_w;
  for(x=0; x<tx; x++,pos++)
      res[pos] = ( 
                  w_side2 * pixels[pos-tx1] + 
                  w_side1 * pixels[pos-tx1] + 
                  w_center* pixels[pos] + 
                  w_side1 * pixels[pos+tx1] + 
                  w_side2 * pixels[pos+tx2] 
                  )/sum_w;
  
  #if defined(USE_OPENMP)
  #pragma omp parallel for
  #endif
  for(y=2; y<ty-2; y++) {
    int x,pos = y*tx;
    for(x=0; x<tx; x++,pos++)
      res[pos] = ( 
                  w_side2 * pixels[pos-tx2] + 
                  w_side1 * pixels[pos-tx1] + 
                  w_center* pixels[pos] + 
                  w_side1 * pixels[pos+tx1] + 
                  w_side2 * pixels[pos+tx2] 
                  )/sum_w;
  }
  pos = (ty-2)*tx;
  for(x=0; x<tx; x++,pos++)
      res[pos] = ( 
                  w_side2 * pixels[pos-tx2] + 
                  w_side1 * pixels[pos-tx1] + 
                  w_center* pixels[pos] + 
                  w_side1 * pixels[pos+tx1] + 
                  w_side2 * pixels[pos+tx1] 
                  )/sum_w;
  for(x=0; x<tx; x++,pos++)
      res[pos] = ( 
                  w_side2 * pixels[pos-tx2] + 
                  w_side1 * pixels[pos-tx1] + 
                  w_center* pixels[pos] + 
                  w_side1 * pixels[pos] + 
                  w_side2 * pixels[pos] 
                  )/sum_w;
}
template<typename type_in_t, typename type_out_t>
void _smooth_7_vert(int tx, int ty, const int w_center, const int w_side1, const int w_side2, const int w_side3, 
                    type_in_t* pixels, type_out_t* res) {
  int x,y,pos=0;
  int sum_w = 2*(w_side1 + w_side2 + w_side3) + w_center;
  if(!sum_w)  sum_w=1;
  const int tx1=tx,tx2=2*tx,tx3=3*tx;
  for(x=0; x<tx; x++,pos++)
      res[pos] = ( 
                  w_side3 * pixels[pos] + 
                  w_side2 * pixels[pos] + 
                  w_side1 * pixels[pos] + 
                  w_center* pixels[pos] + 
                  w_side1 * pixels[pos+tx1] + 
                  w_side2 * pixels[pos+tx2] + 
                  w_side3 * pixels[pos+tx3] 
                  )/sum_w;
  for(x=0; x<tx; x++,pos++)
      res[pos] = ( 
                  w_side3 * pixels[pos-tx1] + 
                  w_side2 * pixels[pos-tx1] + 
                  w_side1 * pixels[pos-tx1] + 
                  w_center* pixels[pos] + 
                  w_side1 * pixels[pos+tx1] + 
                  w_side2 * pixels[pos+tx2] + 
                  w_side3 * pixels[pos+tx3] 
                  )/sum_w;
  for(x=0; x<tx; x++,pos++)
      res[pos] = ( 
                  w_side3 * pixels[pos-tx2] + 
                  w_side2 * pixels[pos-tx2] + 
                  w_side1 * pixels[pos-tx1] + 
                  w_center* pixels[pos] + 
                  w_side1 * pixels[pos+tx1] + 
                  w_side2 * pixels[pos+tx2] + 
                  w_side3 * pixels[pos+tx3] 
                  )/sum_w;
  
  #if defined(USE_OPENMP)
  #pragma omp parallel for
  #endif
  for(y=3; y<ty-3; y++) {
    int x,pos = y*tx;
    for(x=0; x<tx; x++,pos++)
      res[pos] = ( 
                  w_side3 * pixels[pos-tx3] + 
                  w_side2 * pixels[pos-tx2] + 
                  w_side1 * pixels[pos-tx1] + 
                  w_center* pixels[pos] + 
                  w_side1 * pixels[pos+tx1] + 
                  w_side2 * pixels[pos+tx2] + 
                  w_side3 * pixels[pos+tx3] 
                  )/sum_w;
  }
  pos = (ty-3)*tx;
  for(x=0; x<tx; x++,pos++)
      res[pos] = ( 
                  w_side3 * pixels[pos-tx3] + 
                  w_side2 * pixels[pos-tx2] + 
                  w_side1 * pixels[pos-tx1] + 
                  w_center* pixels[pos] + 
                  w_side1 * pixels[pos+tx1] + 
                  w_side2 * pixels[pos+tx2] + 
                  w_side3 * pixels[pos+tx2] 
                  )/sum_w;
  for(x=0; x<tx; x++,pos++)
      res[pos] = ( 
                  w_side3 * pixels[pos-tx3] + 
                  w_side2 * pixels[pos-tx2] + 
                  w_side1 * pixels[pos-tx1] + 
                  w_center* pixels[pos] + 
                  w_side1 * pixels[pos+tx1] + 
                  w_side2 * pixels[pos+tx1] + 
                  w_side3 * pixels[pos+tx1] 
                  )/sum_w;
  for(x=0; x<tx; x++,pos++)
      res[pos] = ( 
                  w_side3 * pixels[pos-tx3] + 
                  w_side2 * pixels[pos-tx2] + 
                  w_side1 * pixels[pos-tx1] + 
                  w_center* pixels[pos] + 
                  w_side1 * pixels[pos] + 
                  w_side2 * pixels[pos] + 
                  w_side3 * pixels[pos] 
                  )/sum_w;
}


/* Smooth an image using a Gaussian filter.
   in-place is ok (img can be ==res)
*/
template<typename type_t>
void _smooth_gaussian_alltype( const int tx, const int ty, type_t* img, float _sigma, type_t* res ) {
  const float MAX_SIGMA = 1.86f;
  
  type_t* img2 = img;
  if(_sigma>MAX_SIGMA) {  // reallocate if more than one smoothing pass is required
    img2 = NEWA(type_t,tx*ty);
    memcpy(img2,img,tx*ty*sizeof(type_t));
  }
  type_t* tmp = NEWA(type_t,tx*ty);
  type_t* old_res = res;
  
  float remaining = _sigma*_sigma;
  while( 1 ) {
    float sigma = MIN(MAX_SIGMA,sqrt(remaining));
    remaining -= sigma*sigma;
    
    // compute gaussian filter coefficients
    const int wcenter = 1000;
    const int wside1 = int(0.5 + wcenter*exp( -pow2(1./sigma)/2 ));
    const int wside2 = int(0.5 + wcenter*exp( -pow2(2./sigma)/2 ));
    const int wside3 = int(0.5 + wcenter*exp( -pow2(3./sigma)/2 ));
    const int wside4 = int(0.5 + wcenter*exp( -pow2(4./sigma)/2 ));
    assert( wside4 < wcenter/10 || !"error: smoothing is too large" );
    
    if ( wside2 < wcenter/10 ) {
      _smooth_3_horiz( tx, ty, wcenter, wside1, img2, tmp );
      _smooth_3_vert(  tx, ty, wcenter, wside1, tmp, res );
    } else if( wside3 < wcenter/10 ) {
      _smooth_5_horiz( tx, ty, wcenter, wside1, wside2, img2, tmp );
      _smooth_5_vert(  tx, ty, wcenter, wside1, wside2, tmp, res );
    } else {
      _smooth_7_horiz( tx, ty, wcenter, wside1, wside2, wside3, img2, tmp );
      _smooth_7_vert(  tx, ty, wcenter, wside1, wside2, wside3, tmp, res );
    }
    
    if(remaining < 0.001)
      break;
    else {
      type_t* tmp3;
      tmp3 = img2;
      img2 = res;
      res = tmp3;
    }
  }
  
  if(res!=old_res) { // copy to true res
    memcpy(old_res,res,tx*ty*sizeof(type_t));
    img2 = res;
  }
  if(_sigma>MAX_SIGMA) 
    free(img2);
  free(tmp);
}


template<typename type_t>
void smooth_gaussian_image(type_t *out, const int h, const int w, const int channels, type_t *in, const float sigma){
    for( int l = 0 ; l < channels ; l++ ){
        _smooth_gaussian_alltype(w, h, in + l*w*h, sigma, out + l*w*h);
    }
}



void rgb_to_lab_inplace(float *image, const int npix, const float color_attenuation){
    float *red = image, *green = image+npix, *blue = image+2*npix;
    const float T=0.008856;
    for(int i = 0 ; i<npix ; ++i){
        float rpix = red[i]/255.0f, gpix = green[i]/255.0f, bpix = blue[i]/255.0f;
        float X=0.412453 * rpix + 0.357580 * gpix + 0.180423 * bpix;
        float Y=0.212671 * rpix + 0.715160 * gpix + 0.072169 * bpix;
        float Z=0.019334 * rpix + 0.119193 * gpix + 0.950227 * bpix;        
        X/=0.950456;
        Z/=1.088754;
        float Y3 = pow(Y,(float) 1./3);      
        float fX = X>T ? pow(X, (float) 1./3) : 7.787 * X + 16/116.;
        float fY = Y>T ? Y3 : 7.787 * Y + 16/116.;
        float fZ = Z>T ? pow(Z, (float) 1./3) : 7.787 * Z + 16/116.;

        red[i] = Y>T ? 116 * Y3 - 16.0 : 903.3 * Y; // L
        green[i] = 500 * (fX - fY); // a
        blue[i] = 200 * (fY - fZ); // b

        if( color_attenuation>0 ) {
            // correct La*b*: dark area or light area have less reliable colors
            float correct_lab = exp(-color_attenuation*pow2(pow2(red[i]/100) - 0.6));
            green[i] *= correct_lab;
            blue[i] *= correct_lab;
        } 
        
    }
}


// hog_params = dict(img_from_lab = True,grad_mode = 0, hog_sharpen = 0, presmooth_sigma = 1.0, mid_smoothing = 1.5, hog_sigmoid = 0.2,  hog_sig_offset = 0, post_smoothing = 1.0, ninth_dim = 0.3, norm_hog = True )
float* extract_hog(const int h, const int w, float* lab){ 
    // img = lab[0]*2.55
    
    const int npix = h*w;
    float *hog = NEWA(float, 9*npix); // output
    memset(hog, 0, 9*npix*sizeof(float));
    uint8_t *L = NEWA(uint8_t, npix); // luminance
    // take luminance
    for(int p=0 ; p<npix ; p++)
        L[p] = (uint8_t) (2.55*lab[p]);
    // presmooth image
    smooth_gaussian_image(L, h, w, 1, L, 1.0); 
    // horizontal gradient
    float *gradx = NEWA(float, npix), *grady = NEWA(float, npix);
    for(int j=0 ; j<h ; j++){
        int p = j*w;
        gradx[p] = L[p+1]-L[p];
        p+=1;
        for(int i=1 ; i<w-1 ; i++,p++){
            gradx[p] = L[p+1] - L[p-1];
        }
        p = j*w+w-1;
        gradx[p] = L[p]-L[p-1];
    }
    // vertical gradient
    { // line 0
        for(int p=0 ; p<w ; p++){
            grady[p] = L[p+w]-L[p];
        }
    }
    for( int j=1 ; j<h-1 ; j++){
        for(int p = j*w ; p<(j+1)*w; p++){
            grady[p] = L[p+w] - L[p-w];
        }
    }
    { // last line
        for(int p=npix-w ; p<npix ; p++){
            grady[p] = L[p]-L[p-w];
        }    
    }
    free(L);
    // build hog using fast cos projection
    const int n_ori = 8;
    for( int l=0 ; l<n_ori ; l++){
        float angle = -2.0f * (l-2) * M_PI / n_ori;
        float kos = cos(angle);
        float zin = sin(angle);
        float *res = hog + l*npix;
        for(int p=0 ; p<npix ; p++){
            float value = kos * gradx[p] + zin * grady[p];
            res[p] = value > 0 ? value : 0;
        }
    }
    free(gradx);
    free(grady);
    // smoothing
    smooth_gaussian_image(hog, h, w, 8, hog, 1.5f); 
    // sigmoid
    for(int p=0 ; p<8*npix ; p++){
        hog[p] = 2.0f / ( 1.0f + expf(-0.2*hog[p])) - 1.0f;
    }    
    // smoothing
    smooth_gaussian_image(hog, h, w, 8, hog, 1.0f); 
    // add ninth dimension and normalize per-pixel
    for( int p=0 ; p<npix ; p++){
        //hog[8*npix+p] = 0.3f;
        float norm = 0.09f;
        for( int l=0 ; l<8 ; l++){
            norm += pow2(hog[l*npix+p]);
        }
        norm = sqrtf(norm)+1e-8;
        for( int l=0 ; l<8 ; l++){
            hog[l*npix+p] /= norm;
        }
        hog[8*npix+p] = 0.3f/norm;
    }    
    return hog;
}


void compute_diff_map(float *out, const int h, const int w, const int channels, float *img1, float *img2, float *flow, const float lr_0, const float lr_1){
    const int npix = w*h;
    bool *invalid = NEWA(bool, npix);
    for(int j=0 ; j<h ; j++){
        for (int i=0 ; i<w ; i++){
            const int p = w*j+i;
            const float xx = i + flow[p], yy = j + flow[p+npix];
            const float x = floor(xx), y = floor(yy);
            invalid[p] = ( xx<0 || xx>w-1.001f || yy<0 || yy>h-1.001f);
            const float dx = xx-x, dy = yy-y;
            const int x1 = MIN(w-1, MAX(0, x)),
                      x2 = MIN(w-1, MAX(0, x+1)),
                      y1 = MIN(h-1, MAX(0, y)),
                      y2 = MIN(h-1, MAX(0, y+1));
            out[p] = 0.0f;
            for( int c=0 ; c<channels ; c++){
                float v1 = img1[p+c*npix], v2 =
                     img2[w*y1+x1+c*npix]*(1.0f-dx)*(1.0f-dy) + 
                     img2[w*y1+x2+c*npix]*dx*(1.0f-dy) + 
                     img2[w*y2+x1+c*npix]*(1.0f-dx)*dy + 
                     img2[w*y2+x2+c*npix]*dx*dy;
                out[p] += pow2(v1-v2);
            }
            out[p] = sqrt(out[p]);
        }
    }  
    smooth_gaussian_image(out, h, w, 1, out, 1.0);
    for( int p=0 ; p<npix ; p++){
        out[p] = invalid[p] ? 0.0f : 1.0f/(1.0f+exp(-lr_0-lr_1*out[p]));
    }
    free(invalid);
}

void image_warping_error(float *out_lab, float *out_hog, const int h, const int w, float *image1, float *image2, float *flow){

    // convert to Lab
    rgb_to_lab_inplace(image1, h*w, 10.0);
    rgb_to_lab_inplace(image2, h*w, 10.0);   
    
    // extract hog
    float *hog1 = extract_hog(h, w, image1);
    float *hog2 = extract_hog(h, w, image2);
    
    // smooth images
    smooth_gaussian_image(image1, h, w, 3, image1, 1.0);
    smooth_gaussian_image(image2, h, w, 3, image2, 1.0);
        
    // compute warping
    compute_diff_map(out_lab, h, w, 3, image1, image2, flow, -2.0f, 0.15f);
    compute_diff_map(out_hog, h, w, 9, hog1, hog2, flow, -3.5f, 5.0f);

    free(hog1);
    free(hog2);
}

float *input3darray_to_c_order_layer(const mxArray *p){
    const int *dims = mxGetDimensions(p);
    const int h = dims[0], w = dims[1], channels = dims[2];
    float *out = NEWA(float, h*w*channels);
    float *in = (float*) mxGetData(p);
    for( int c = 0 ; c < channels ; c++ ) {
        float *inptr = in + c * w * h;
        float *outptr = out + c * w * h;
        for ( int j = 0 ; j<h ; j++ ){
            for ( int i = 0 ; i<w ; i++ ){
                outptr[j*w+i] = inptr[i*h+j];
            }
        }
    }
    return out;
}

void ptr_to_outputarray(mxArray *p, float *data){
    const int h = mxGetM(p), w = mxGetN(p);
    float *array_data = (float*) mxGetData(p);
    for( int j=0 ; j<h ; j++ ){
        for( int i=0 ; i<w ; i++){
            array_data[i*h+j] = data[j*w+i];
        }
    }
}

void mexFunction( int nl, mxArray *pl[], int nr, const mxArray *pr[] ) {
    
    assert ( nl == 2);
    assert ( nr == 3); 
    
    // The code is originally written for C-order arrays.
    // We thus transpose all arrays in this mex-function which is not efficient...
    
    const int *pDims;
    assert( mxGetNumberOfDimensions(pr[0]) == 3);
    assert( mxIsClass(pr[0], mxSINGLE_CLASS));
    pDims = mxGetDimensions(pr[0]);
    assert( pDims[2]==3);
    const int h = pDims[0], w = pDims[1];
    float *image1 = input3darray_to_c_order_layer( pr[0] );
   
    assert( mxGetNumberOfDimensions(pr[1]) == 3);
    assert( mxIsClass(pr[1], mxSINGLE_CLASS));
    pDims = mxGetDimensions(pr[1]);
    assert( pDims[0]==h && pDims[1]==w && pDims[2]==3 );
    float *image2 = input3darray_to_c_order_layer( pr[1] );

    assert( mxGetNumberOfDimensions(pr[2]) == 3);
    assert( mxIsClass(pr[2], mxSINGLE_CLASS));
    pDims = mxGetDimensions(pr[2]);
    assert( pDims[0]==h && pDims[1]==w && pDims[2]==2 );
    float *flow = input3darray_to_c_order_layer( pr[2] );


    float *out_lab = NEWA(float, h*w);
    float *out_hog = NEWA(float, h*w);

    
    image_warping_error(out_lab, out_hog, h, w, image1, image2, flow);


    pl[0] = mxCreateNumericMatrix(h, w, mxSINGLE_CLASS, mxREAL);
    pl[1] = mxCreateNumericMatrix(h, w, mxSINGLE_CLASS, mxREAL);
    ptr_to_outputarray(pl[0], out_lab);
    ptr_to_outputarray(pl[1], out_hog);

    free(image1);
    free(image2);
    free(flow);
    free(out_lab);
    free(out_hog);
    

}
