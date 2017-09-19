#include <stdio.h>
#include <math.h>
#include <string.h>
#define forall(range_h) for(i=0; i<range_h; i++)
#define copy_vector(y,x,len) memcpy(y,x,sizeof(double)*len)

void scale_vector(double *y, const double alpha, const double *x, const int len)
{
	for(int i=0; i<len; i++) y[i] = alpha*x[i];
}


void quick_sort (double *a, int n) {
    int i, j;
    double p, t;
    if (n < 2)
        return;
    p = a[n / 2];
    for (i = 0, j = n - 1;; i++, j--) {
        while (a[i] < p)
            i++;
        while (p < a[j])
            j--;
        if (i >= j)
            break;
        t = a[i];
        a[i] = a[j];
        a[j] = t;
    }
    quick_sort(a, i);
    quick_sort(a + i, n - i);
}

 double norm_one(double *x, int len)
{
  double n = 0.0;
  for(int i=0; i<len; i++) n = x[i] > 0 ? n+x[i] : n-x[i];
    return n;
}
 
 
void proj_normball_one(double *xproj, double *x, const double c, const int len)
{
    double tmp_norm, tmp_lambda, tmp_sum, tmp_glambda;
    int i;
    //int len1=len+1;
    
    
    scale_vector(xproj, 1/c, x, len);
    printf("xproj[0] is %f and xproj[1] is %f \n", xproj[0], xproj[1]);
    tmp_norm=norm_one(xproj,len);
    
    printf("tmp_norm is %f \n", tmp_norm);
    
    if(tmp_norm<1)
    {
        tmp_lambda=0;
        printf("it sets lambda equal to zero. \n");
    }
    else
    {
        
        
        /* Next we find |a_k| */
        forall(len) xproj[i] = (xproj[i] < 0.0) ? -xproj[i] : xproj[i];
        
        /*We sort in ascending order */
        quick_sort (xproj, len);
        printf("After sorting xproj[0] is %f and xproj[1] is %f \n", xproj[0], xproj[1]);
        /* Compute the table from bottom of Algorithm I */
        
        tmp_sum=-1;
        for(i=len-1; i>0; i--){
         
            tmp_sum=tmp_sum+xproj[i];
            
            printf("The first i is %d \n",i);
            
            printf("temporary sum value is %f \n", tmp_sum);
            tmp_glambda=(i-len)*xproj[i-1]+tmp_sum;
            
            printf("tmp_glambda is %f \n", tmp_glambda);
            if (tmp_glambda>0){
                
                break;                
            }
            
        }
        
        printf("The loop is broken at the iteration %d \n",i);
        /* if tmp_lambda >0 then oir k=i else our k=0 i.e. the sign will change for k=0 */
        if(tmp_glambda>0)
        
        {
            printf("it goes with tmp_lambda is greater than zero \n");
            tmp_lambda=tmp_sum/(len-i);
        
        }
        
        else
        {   
            
            printf("it goes i equal to zero iteration \n");
            tmp_lambda=(tmp_sum+xproj[0])/len;
            
        }
        
        printf("The value of the lambda is %f \n", tmp_lambda);
        /* Compute the projection using equation 8. Here I am using macro inside macro */
        forall(len) xproj[i] = (x[i]+tmp_lambda <= 0.0) ? x[i]+tmp_lambda : ((x[i]-tmp_lambda >= 0)? x[i]-tmp_lambda: 0);
    }
}
        
    
    


int main () {
    double a[] = {4.5, 65.2, 65.25, -31, 0, 99.2, 2, 83, 782, 1};
    int n = sizeof a / sizeof a[0];
    int i;
    
    // Sorting test
    
    
    for (i = 0; i < n; i++)
        printf("%f%s", a[i], i == n - 1 ? "\n" : " ");
    quick_sort(a, n);
    for (i = 0; i < n; i++)
        printf("%f%s", a[i], i == n - 1 ? "\n" : " ");
        
        
        for(i=n-1; i>0; i--){
         
            printf("%d \n",i);
            
        }
        
        printf("%d \n",i);
        
        
        // double macro test
        double b[]={-1, 2, 5};
        
        
        forall(3) b[i] = (b[i] < 0.0) ? -b[i] : ((b[i]<3)? b[i]: -b[i]);
        
        printf("%f and %f and %f \n",b[0], b[1], b[2]);
        
        
        // Normball test
        
        double c[]={-7, -5, -6};
        const int len=3;
        const double c_tmp=1;
        
        double d[3];
        
        proj_normball_one(d, c,  c_tmp, len);
        
        printf("d1 is %f and d2 is %f and d3 is %f \n",d[0],d[1], d[2]);
    return 0;
}