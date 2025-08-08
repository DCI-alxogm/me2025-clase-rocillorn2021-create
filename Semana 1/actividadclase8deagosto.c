//8 de Agosto del 2025
/*Hacer un programa que imprima 5 nùmeros flotantes, y el valor de 
las siguientes funciones. */

/* f(x)= x2
/ f(x)= log (x+1)
/f (x) = exp (x) + cos (x) */

/*Imprime 5 nùmeros flotantes*/

#include <stdio.h>
#include <math.h>

int main (){
  float valores[5];
  valores [0] = 1.5;
  valores [1] = 2.5;
  valores [2] = 3.5;
  valores [3] = 4.5;
  valores [4] = 5.5;

  printf("Impresiòn 5 nùmeros floltante.");
printf("Impresiòn 5 nùmeros floltante son: %f %f %f %f %f\n",  valores [0],  valores [1], valores [2], valores [3], valores [4]);
  


  /*Imprime la funciones 

  / f(x)= x2
  / f(x)= log (x+1)
  /f (x) = exp (x) + cos (x) */

  float x;
  x=6.5;

  printf("f(x) = x^2 = %f \n" , x*x);
  printf("f(x) = log(x+1 = %f \n" , log(x+1));
  printf("f(x) = (exp(x) + cos(x)) = %f \n" , (exp(x) + cos(x)));
  }