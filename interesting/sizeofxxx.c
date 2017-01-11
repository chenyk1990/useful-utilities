#include<stdio.h>
/*the following output will be system dependent*/
/*gcc sizeofxxx.c && ./a.out*/

int main()
{
 int i;
 printf("size of (int) is %lu \n",sizeof(int)); /*lu: unsigned long*/
 printf("size of (int) is %d \n",sizeof(int));  /*with warning*/
 printf("size of (int) is %lu \n",sizeof(long int)); 
 printf("size of (int) is %lu \n",sizeof(unsigned int));  
 printf("size of (int) is %lu \n",sizeof(unsigned long));
 printf("size of (int) is %lu \n",sizeof(unsigned));
 printf("size of (int) is %lu \n",sizeof(long));
 printf("size of (int) is %lu \n",sizeof(float));
 printf("size of (int) is %lu \n",sizeof(double));
 printf("size of (int) is %lu \n",sizeof(long double));    
 printf("size of (int) is %lu \n",sizeof(char));    
 return 0;
}