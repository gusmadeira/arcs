/////Esse programa faz a conversao de coordenada cartesianas dadas em relação a 
/////a um corpo central  achatado para elementos orbitais geometricos.

//Programador: André Izidoro Ferreira da Costa


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#define nlc 4

//nlc Número de linhas do cabeçalho
char *subs_str(char *string, char *substring, char *nova)
  {
    char *extrai;
    int tamanho1,tamanho2,contador;

    tamanho1 = strlen(substring);
    tamanho2 = strlen(nova);

    if((tamanho1 > tamanho2) || (tamanho2 > tamanho1))
      return(" ");
    else
      {
        extrai = strstr(string,substring);

        if(extrai)
          {
            for(contador = 0;contador < tamanho1; contador++)
              string[(extrai - string) + contador] = nova[contador];
            return(string);
          }
        else
         return(" ");
       }
  }

//ONDE E COMO ATRIBUO O NOME DO ARQUIVO DE ENTRADA DA VARIÁVEL ARGV[1]?
int main( int argc, char *argv[])
{
	
	FILE *ArqE;
	FILE *ArqS;

	char Nome[20];
	char Linha[200];

	
        ArqE = fopen(argv[1], "r");

	strcpy(Nome,argv[1]);
//	perror("Passou por aqui");
	strcat(Nome,"_out");





	ArqS = fopen(Nome,"w");

	

	 double sem,exc,INC,inclina;
	 double x,y,z,vx,vy,vz,MU;
	 double J2,J4,J6,RS,MS,mass;
	 double Gauss2,SMA,ECC,wpi;
	 double NMEDI,w,id,modulo;
	 double kappa,eta2,chi2,nu;
	 double alpha,alpha2,longnodo;
	 double PI,lambda,epson,ai;
	 double L,Ldot,zdot,M,lom,R;
	 double Rc, Lc, zc, Rcdot, Lcdot,zcdot;
	 double Erre,aux;

	 double Erredot,Pos1,Pos2,Pos3,Vel1,Vel2,Vel3,an,en;
	 double alphaa[3];
	 int itera,i;

//	char  ARQCI[20];
	PI = 3.1415926535897932384626433832795;
	Gauss2=6.6743e-11;
	J2=16290.543820e-6;
	J4=-936.700366e-6;
	J6=0.0E0;
	RS=60330e3;
	MS=5.68683765495e26;

//	fscanf(ArqE,"%lE %lE %lE %lE %lE %lE %lE %lE %lE",&id,&R,&Pos1,&Pos2,&Pos3,&Vel1,&Vel2,&Vel3,&mass);	
	fscanf(ArqE,"%lE %lE %lE %lE %lE %lE %lE %lE %lE %lE",&id,&Pos1,&Pos2,&Pos3,&Vel1,&Vel2,&Vel3,&mass,&an,&en);	

	do{


	
	Erre=pow((Pos1*Pos1 + Pos2*Pos2),0.5);

	z=Pos3;
	zdot=Vel3;
	
	if(Pos1< 0.0)
	L=atan(Pos2/Pos1)+1.0*PI;
        else if (atan(Pos2/Pos1)< 0.0)
	L=atan(Pos2/Pos1)+2.0*PI;
	else if ((Pos1==0.0) && (Pos2>0.0))
	L=(PI/2.0);
	else if ((Pos1==0.0) && (Pos2<0.0))
	L=((3.0*PI)/2.0);
	else
	L=atan(Pos2/Pos1);	
	


	do{ 
	if (L>2.0*PI)
	L=L-2*PI;
	else if (L < 0.0)
	L=L+2*PI;
	}
	while(L > (2.0*PI)|| L< 0.0);



	Erredot=Vel1*cos(L)+Vel2*sin(L);
	Ldot=(-Vel1*sin(L)+Vel2*cos(L))/Erre;



	itera=0;



 //! initialization of variables
      
	SMA=Erre;
        ECC=0.0;
        INC=0.0;
	
	Rc=0.0; 
	Lc=0.0;
	zc=0.0;
	Rcdot=0.0;
	Lcdot=0.0;
	zcdot=0.0;
		
	






        MU=Gauss2*(MS+mass);


	 epson=1.0E-15;

       do{

	ai=SMA;

	NMEDI=(pow(MU/pow(SMA,3),(1.0/2.0)))*(1+((3.0/4.0)*J2*
	pow((RS/SMA),2))-((15.0/16.0)*J4*pow((RS/SMA),4))+((35.0/32.0)*J6*
	pow((RS/SMA),6))-((9.0/32.0)*pow(J2,2)*pow((RS/SMA),4))+((45.0/64.0)*J2*J4*
	pow((RS/SMA),6))+((27.0/128.0)*pow(J2,3)*pow((RS/SMA),6))+(3*J2*
	pow(ECC,2)*pow((RS/SMA),2))-(12*pow((RS/SMA),2)*J2*pow(INC,2)));



//Clculo da freq?ncia epicclica horizontal
	kappa=(pow(MU/pow(SMA,3),(1.0/2.0)))*(1-((3.0/4.0)*J2*
	pow(RS/SMA,2))+
	((45.0/16.0)*J4*pow((RS/SMA),4))-((175.0/32.0)*J6*pow(RS/SMA,6))-
	((9.0/32.0)*
	pow(J2,2)*pow((RS/SMA),4))+((135.0/64.0)*J2*J4*pow(RS/SMA,6))-
	((27.0/128.0)
	*pow(J2,3)*pow(RS/SMA,6))-(9*pow(RS/SMA,2)*J2*pow(INC,2)));    



//Clculo da freq?ncia epclica vertical




	nu=(pow(MU/pow(SMA,3),(1.0/2.0)))*(1+((9.0/4.0)*J2*pow(RS/SMA,2))-
	((75.0/16.0)*J4*pow((RS/SMA),4))+((245.0/32.0)*J6*pow(RS/SMA,6))-
	((81.0/32.0)
	*pow(J2,2)*pow((RS/SMA),4))+((675.0/64.0)*J2*J4*pow(RS/SMA,6))+
	((729.0/128.0)*pow(J2,3)*pow(RS/SMA,6))+
	6*pow(RS/SMA,2)*J2*pow(ECC,2)-
	(51.0/4.0)*pow(RS/SMA,2)*J2*pow(INC,2));  



	eta2=((MU/pow(SMA,3)))*(1-(2*J2*pow(RS/SMA,2))+
	((75.0/8.0)*J4*pow((RS/SMA),4))-((175.0/8.0)*J6*pow(RS/SMA,6)));


	chi2=(MU/pow(SMA,3))*(1+(15.0/2.0)*J2*pow(RS/SMA,2)-
	(175.0/8.0)*pow((RS/SMA),4)*J4+(735.0/16.0)*pow(RS/SMA,6)*J6);




      



	alphaa[1]=1.0/3.0*(2*nu+kappa);
	alphaa[2]=2*nu-kappa;
	alpha2=alphaa[1]*alphaa[2];
	alpha=sqrt(alpha2);



	sem=(Erre-Rc)/(1-(Ldot-Lcdot-NMEDI)/(2*NMEDI));


	exc=sqrt(pow(((Ldot-Lcdot-NMEDI)/(2*NMEDI)),2)+
	pow(((Erredot-Rcdot)/(sem*kappa)),2));
	

	inclina=pow((pow(((z-zc)/(sem)),2)+pow(((zdot-zcdot)/(sem*nu)),2)),(1.0/2.0));


	do{ 
	if(inclina>PI)
	inclina=inclina-PI;
	else if (inclina < 0.0)
	inclina=inclina+PI;
	}while(L< 0.0 || inclina>PI);


	lambda=L-Lc-2*(NMEDI/kappa)*((Erredot-Rcdot)/(sem*kappa));
	
	

	do{ 
	if(lambda> 2*PI) 
	lambda=lambda-2*PI;
	else if(lambda < 0.0)
	lambda=lambda+2*PI;
	}
	while((lambda> (2*PI))||lambda <0.0);
	
	 

	


	M=atan((Erredot-Rcdot)/(sem*kappa*(1-((Erre-Rc)/sem))));

	if((-Erre+Rc+sem)/(sem*exc)<0.0)
	M=M+1.0*PI;
	else if(M<0.0)
	M=M+2*PI;
	else if (((-Erre+Rc+sem)/(sem*exc))==00 && (Erredot-Rcdot)/(sem*kappa*kappa)>0.0)
	M=(PI/2.0);
	else if (((-Erre+Rc+sem)/(sem*exc))==00 && (Erredot-Rcdot)/(sem*kappa*kappa)<0.0)
	M=((3.0*PI)/2.0);
		
		
	do{ 
	if(M> 2*PI) 
	M=M-2*PI;
	else if(M< 0.0)
	M=M+2*PI;
	}
	while((M> (2*PI))||M <0.0);
	
	
	
//	    lom=lambda; // sabemos que Omega= 0  criterio escolhido quando geramos os dados (caso plano)
                  // lom é definido como lom=lamda - Omega
			       // e lambda como lambda= M + w 
			       //Para mais veja (Renner&Sicardy 2006, Celestial Mech. Dyn. Astr.)
			    



	lom=atan((nu*(z-zc))/(zdot-zcdot));
	


	if(((zdot-zcdot)/(sem*inclina*nu) <0.0))
	lom=lom+PI;	
	else if(lom<0.0)
	lom=lom+2*PI;
	else if((((nu*(z-zc))/(sem*inclina))>0.0 && ((zdot-zcdot)/(sem*inclina*nu)==0.0)))
	lom=PI/2.0;	
	else if((((nu*(z-zc))/(sem*inclina))<0.0 && ((zdot-zcdot)/(sem*inclina*nu)==0.0)))
	lom=(3.0*PI)/2.0;


  	do{ 
	if(lom> 2*PI) 
	lom=lom-2*PI;
	else if(lom< 0.0)
	lom=lom+2*PI;
	}
	while((lom>(2*PI))||lom <0.0);

	


	SMA=sem;
	ECC=exc;
	INC=inclina;

	Rc=SMA*pow(ECC,2)*((3.0/2.0)*(eta2/pow(kappa,2.0))-1-(1.0/2.0)
	*((eta2/(pow(kappa,2.0)))
	*cos(2*M)))+SMA*pow(INC,2)*((3.0/4.0)*(chi2/pow(kappa,2))-1+(chi2/
	(4*alpha2)*cos(2*lom)));

	Lc=(NMEDI/kappa)*pow(ECC,2)*((3.0/4.0)+eta2/(2*pow(kappa,2)))*sin(2*M)
	-pow(INC,2)*(chi2/(4*alpha2))*(NMEDI/nu)*sin(2*lom);

	zc=SMA*INC*ECC*((chi2/(2*kappa*alphaa[1]))*sin(lom+M)-(3.0/2.0)
	*((chi2)/(kappa*alphaa[2]))*sin(lom-M));

	Rcdot=SMA*pow(ECC,2)*(eta2/kappa)*sin(2*M)-SMA*pow(INC,2)*
	((chi2/(2*alpha2)))*nu*sin(2*lom);
     

	Lcdot=NMEDI*pow(ECC,2)*((7.0/2.0)-3*(eta2/pow(kappa,2))-(pow(kappa,2)
	/(2*pow(NMEDI,2)))+((3.0/2.0)+(eta2/pow(kappa,2)))*cos(2*M))+pow(INC,2)*
	NMEDI*(2-pow(kappa,2)/(2*pow(NMEDI,2))-(3.0/2.0)*(chi2/pow(kappa,2)
	)-(chi2/(2*alpha2))*cos(2*lom));
	
	zcdot=SMA*INC*ECC*(((chi2*(kappa+nu))/(2*kappa*alphaa[1]))
	*cos(lom+M)+(3.0/2.0)*((chi2*(kappa-nu))/
	(kappa*alphaa[2]))*cos(lom-M));



	modulo=(SMA-ai);

	if(modulo<0.0)
	modulo=modulo*(-1.0);

	if ((modulo<epson||itera==999) && !feof(ArqE)){
		
	longnodo=lambda-lom;

	wpi=lambda-M;
			
  		  	do{ 
			if(longnodo>2*PI) 
			longnodo=longnodo-2*PI;
			else if(longnodo< 0.0)
			longnodo=longnodo+2*PI;
			}while((longnodo> (2*PI))||longnodo < 0.0);


			do{ 
			if(wpi> 2*PI) 
			wpi=wpi-2*PI;
			else if(wpi< 0.0)
			wpi=wpi+2*PI;
			}while((wpi> (2*PI))||wpi < 0.0);
	
	
	fprintf(ArqS,"% 10.5lf % 18.15lf % 18.15lf % 18.15lf  %18.15lf  %18.15lf  %18.15lf  %10.5lf\n",id,sem,exc,inclina,wpi*(180.0/PI),longnodo*(180.0/PI),lambda*(180.0/PI),mass);
//	fprintf(ArqS,"  %018.15lf  %018.15lf  %018.15lf  %019.15lf\n",id,NMEDI-nu,NMEDI,mass);
//	fprintf(ArqS,"  %019.15lf    %019.15lf   %019.15lf   %019.15lf    %019.15lf   %019.15lf   %019.15lf   %019.15lf \n",id,sem,exc,inclina*(180.0/PI),longnodo*(180.0/PI),wpi*(180.0/PI),lambda*(180.0/PI),mass);
	}
	

	

//	printf("itera= %d, modulo=%le\n",itera,modulo);

	itera++;

	
	}
	while (modulo>epson&&itera<1000);

//	fscanf(ArqE,"%lE %lE %lE %lE %lE %lE %lE %lE %lE",&id,&R,&Pos1,&Pos2,&Pos3,&Vel1,&Vel2,&Vel3,&mass);
	fscanf(ArqE,"%lE %lE %lE %lE %lE %lE %lE %lE  %lE %lE",&id,&Pos1,&Pos2,&Pos3,&Vel1,&Vel2,&Vel3,&mass,&an,&en);
}
while(!feof(ArqE));
		fclose(ArqE);
		fclose(ArqS);


return 0;
}
