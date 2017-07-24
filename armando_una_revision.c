#include <stdio.h>
#include <stdlib.h>
#include <math.h>



struct Par {double alfa; double gamma1; double gamma2;double p;double p1;double p2;} aa;


void takens(int n,double v[],double dv[],double t) {
                
    double x,y,z,w;
    x=v[0]; y=v[1]; z=v[2];w=v[3];

   
    
 dv[0]=  y;
 dv[1]=  -1*(aa.gamma1*aa.gamma1)*(x+0.8*x*x*x)+3*aa.gamma1*(10*aa.p-x*x)*y; //labio
    
    
 dv[2]=  w;
 dv[3]=  -1*(aa.gamma2*aa.gamma2)*(z+0.8*z*z*z)+aa.p2*aa.gamma2*(10*aa.p1-z*z)*w; //ventanita

        return;
}



int main(){

  int i;
  double t,tiempot,taux;
  double dt ;
  double v[4];
 
  		//condiciones iniciales
  v[0]=.00001; v[1]=0.001;
    v[2]=0.001;v[3]=0.001;
    
    aa.gamma1=5000; aa.gamma2=4000;
    
    tiempot=0;
    dt=1/882000.0;
    i=0;
    taux=0;
    
    while (tiempot<3.)	{
        
        tiempot=tiempot+dt;
      /*
        if((tiempot>1.0)&&(tiempot<1.2)){
            aa.p=0.1+0.25*sin(3.14159*(tiempot-1.0)/.40)-0.25*v[2];
            aa.p1=-0.5+0.35*sin(3.14159*(tiempot/1.0)/.40)+0.1*v[0];
            aa.p2=9.0;
            rk4(takens,v,4,tiempot,dt);
            if(taux==20) {
                printf("%lg\t%lg\t%lg\t%lg\n",tiempot,v[0],v[2],v[0]+0.1*v[2]);
                taux=0;}
            taux++;
        }*/
        
        
        if((tiempot>1.20)&&(tiempot<=1.4)){
            aa.p=0.1+0.25*sin(3.14159*(tiempot-1.2)/1.80)-0.25*v[2];
            aa.p1=-0.1+.25*sin(3.14159*(tiempot-1.2)/1.80)+0.1*v[0];
            aa.p=0.35;
            aa.p1=0.35;
           
            aa.p2=3.50;
            rk4(takens,v,4,tiempot,dt);
            if(taux==20) {
                printf("%lg\t%lg\t%lg\t%lg\t%lg\n",tiempot,v[0],v[2],sin(3.14159*(tiempot-1.2)/0.2)*(v[0]+0.1*v[2]),0.1+0.25*sin(3.14159*(tiempot-1.2)/1.80));
                taux=0;}
            taux++;
        }
        
        if((tiempot>1.4)&&(tiempot<1.6)){
            aa.p=-0.5;
            aa.p1=-0.5;
            aa.p2=3.50;
            rk4(takens,v,4,tiempot,dt);
            if(taux==20) {
                printf("%lg\t%lg\t%lg\t%lg\t%lg\n",tiempot,v[0],v[2],v[0]+0.1*v[2],aa.p);
                taux=0;}
            taux++;
        }
        
        if((tiempot>1.6)&&(tiempot<1.8)){
            aa.p=0.1+0.25*sin(3.14159*(tiempot-1.6)/1.80)-0.25*v[2];
            aa.p1=-0.5;
            aa.p2=3.50;
            rk4(takens,v,4,tiempot,dt);
            if(taux==20) {
                printf("%lg\t%lg\t%lg\t%lg\t%lg\n",tiempot,v[0],v[2],sin(3.14159*(tiempot-1.4)/0.2)*(v[0]+0.1*v[2]),0.1+0.25*sin(3.14159*(tiempot-1.6)/1.80));
                taux=0;}
            taux++;
        }
       
        
       
        
        
    }
    
   
    
    
    

}
