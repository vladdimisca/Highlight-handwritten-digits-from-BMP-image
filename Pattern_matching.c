#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct{
               int x_sus,y_sus,x_jos,y_jos;
              }fereastra;

typedef struct{
               double corelatie;
               fereastra punct;
               int nr_culoare;
              }detectie;

void grayscale_image(char* nume_fisier_sursa,char* nume_fisier_destinatie)
{
   FILE *fin, *fout;
   unsigned int latime_img, inaltime_img;
   unsigned char pRGB[3],aux;

   fin = fopen(nume_fisier_sursa, "rb");
   if(fin == NULL)
   	{
   		printf("nu am gasit imaginea sursa din care citesc");
   		return;
   	}

   fout = fopen(nume_fisier_destinatie, "wb+");
   if(fout == NULL)
   	{
   		printf("nu am gasit fisierul in care copiez imaginea grayscale");
   		return;
   	}

   fseek(fin, 18, SEEK_SET);
   fread(&latime_img, sizeof(unsigned int), 1, fin);
   fread(&inaltime_img, sizeof(unsigned int), 1, fin);

   //copiaza octet cu octet imaginea initiala in cea noua
	fseek(fin,0,SEEK_SET);
	unsigned char c;
	while(fread(&c,1,1,fin)==1)
	{
		fwrite(&c,1,1,fout);
		fflush(fout);
	}
	fclose(fin);

	//calculam padding-ul pentru o linie
	int padding;
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;

	fseek(fout, 54, SEEK_SET);
	int i,j;
	for(i = 0; i < inaltime_img; i++)
	{
		for(j = 0; j < latime_img; j++)
		{
			//citesc culorile pixelului
			fread(pRGB, 3, 1, fout);
			//fac conversia in pixel gri
			aux = 0.299*pRGB[2] + 0.587*pRGB[1] + 0.114*pRGB[0];
			pRGB[0] = pRGB[1] = pRGB[2] = aux;
        	fseek(fout, -3, SEEK_CUR);
        	fwrite(pRGB, 3, 1, fout);
        	fflush(fout);
		}
		fseek(fout,padding,SEEK_CUR);
	}
	fclose(fout);
}

void template_matching(char *nume_imagine,char *nume_sablon,double ps, detectie **D,int *nr_detectii,int nr_sablon)
{
    FILE *fi=fopen(nume_imagine,"rb");

    if(fi==NULL)
    {
        printf("Nu am putut deschide fisierul care contine imaginea %s !\n",nume_imagine);
        return;
    }

    FILE *fs=fopen(nume_sablon,"rb");

    if(fs==NULL)
    {
        printf("Nu am putut deschide fisierul care contine sablonul %s !\n",nume_sablon);
        return;
    }

    unsigned int w_sablon,h_sablon,w_img,h_img,padd_sablon,padding;

    fseek(fs,18,SEEK_SET);

    fread(&w_sablon,sizeof(int),1,fs);
    fread(&h_sablon,sizeof(int),1,fs);

    if(w_sablon%4!=0)
        padd_sablon=4-(3*w_sablon)%4;
    else
        padd_sablon=0;

    fseek(fi,18,SEEK_SET);

    fread(&w_img,sizeof(int),1,fi);
    fread(&h_img,sizeof(int),1,fi);

    if(w_img%4!=0)
        padding=4-(3*w_img)%4;
    else
        padding=0;

    int i,j,k,t;
    unsigned char **I,**S;

    I=(unsigned char **)malloc(h_img*sizeof(char *));

    if(I==NULL)
    {
        printf("Eroare la alocarea matricei care contine imaginea!\n");
        return ;
    }

    for(i=0;i<h_img;i++)
    {
        I[i]=(unsigned char *)malloc(3*w_img);

        if(I[i]==NULL)
        {
            printf("Eroare la alocarea liniei %d din matricea care contine imaginea!\n",i);
            return ;
        }
    }

    fseek(fi,54,SEEK_SET);

    for(i=h_img-1;i>=0;i--)
    {
        fread(I[i],3*w_img,1,fi);
        fseek(fi,padding,SEEK_CUR);
    }

    S=(unsigned char **)malloc(h_sablon*sizeof(char *));

    if(S==NULL)
    {
        printf("Eroare la alocarea matricei care contine sablonul %d!\n",nr_sablon);
        return ;
    }

    for(i=0;i<h_sablon;i++)
    {
        S[i]=(unsigned char *)malloc(3*w_sablon);

        if(S[i]==NULL)
        {
            printf("Eroare la alocarea liniei %d din matricea care contine sablonul %d!\n",i,nr_sablon);
            return ;
        }
    }

    fseek(fs,54,SEEK_SET);

    for(i=h_sablon-1;i>=0;i--)
        {
            fread(S[i],3*w_sablon,1,fs);
            fseek(fs,padd_sablon,SEEK_CUR);
        }

    unsigned int n;
    double corelatie,S_fI,S_sb,S_bar,fI_bar;

    n=w_sablon*h_sablon;

    for(i=0;i<h_img-h_sablon+1;i++)
        for(j=0;j<w_img-w_sablon+1;j++)
        {

            corelatie=0.0;
            fI_bar=0.0;
            S_bar=0.0;
            S_sb=0.0;
            S_fI=0.0;

            for(k=0;k<h_sablon;k++)
                for(t=0;t<w_sablon;t++)
                {
                    fI_bar+=(double)I[i+k][3*(j+t)];
                    S_bar+=(double)S[k][3*t];
                }

            S_bar/=n;
            fI_bar/=n;

            for(k=0;k<h_sablon;k++)
                for(t=0;t<w_sablon;t++)
            {
                S_sb+=((double)S[k][3*t]-S_bar)*((double)S[k][3*t]-S_bar);
                S_fI+=((double)I[i+k][3*(j+t)]-fI_bar)*((double)I[i+k][3*(j+t)]-fI_bar);
            }

            S_sb/=n-1;
            S_sb=sqrt(S_sb);

            S_fI/=n-1;
            S_fI=sqrt(S_fI);

            for(k=0;k<h_sablon;k++)
                for(t=0;t<w_sablon;t++)
                corelatie+=((double)I[i+k][3*(j+t)]-fI_bar)*((double)S[k][3*t]-S_bar);

            corelatie/=n*S_fI*S_sb;

            if(corelatie>ps){
                              *D=(detectie *)realloc(*D,(*nr_detectii+1)*sizeof(detectie));

                              if(*D==NULL)
                              {
                                  printf("Eroare la realocarea memoriei pentru vectorul de detectii!\n");
                                  return;
                              }

                              (*D)[*nr_detectii].corelatie=corelatie;
                              (*D)[*nr_detectii].punct.x_sus=i;
                              (*D)[*nr_detectii].punct.y_sus=j;
                              (*D)[*nr_detectii].punct.x_jos=i+h_sablon-1;
                              (*D)[*nr_detectii].punct.y_jos=j+w_sablon-1;
                              (*D)[*nr_detectii].nr_culoare=nr_sablon;

                              (*nr_detectii)++;
                            }
        }

    fclose(fi);
    fclose(fs);

    for(i=0;i<h_img;i++)
        free(I[i]);
    free(I);

    for(i=0;i<h_sablon;i++)
        free(S[i]);
    free(S);

}

void colorare(char *nume_imagine,fereastra f,unsigned char *C)
{
    FILE *fin=fopen(nume_imagine,"rb+");

    if(fin==NULL)
    {
        printf("Eroare la deschiderea fisierului care contine imaginea color !\n");
        return;
    }

    unsigned int w_img,h_img,padding;
    int i;

    fseek(fin,18,SEEK_SET);

    fread(&w_img,sizeof(int),1,fin);
    fread(&h_img,sizeof(int),1,fin);

    if(w_img%4!=0)
        padding=4-(3*w_img)%4;
    else
        padding=0;

    unsigned char **I;

    I=(unsigned char **)malloc(h_img*sizeof(char *));

    if(I==NULL)
    {
        printf("Eroare la alocarea memoriei pentru matricea care retine imaginea in cadrul colorarii!\n");
        return;
    }

    for(i=0;i<h_img;i++)
    {
        I[i]=(unsigned char *)malloc(3*w_img);
        if(I[i]==NULL)
        {
            printf("Eroare la alocarea memoriei pentru linia %d din matricea care retine imaginea!\n",i);
            return;
        }
    }

    fseek(fin,54,SEEK_SET);

    for(i=h_img-1;i>=0;i--)
    {
        fread(I[i],3*w_img,1,fin);
        fseek(fin,padding,SEEK_CUR);
    }

    for(i=f.x_sus;i<=f.x_jos;i++)
    {
        I[i][3*f.y_sus]=I[i][3*f.y_jos]=C[2];
        I[i][3*f.y_sus+1]=I[i][3*f.y_jos+1]=C[1];
        I[i][3*f.y_sus+2]=I[i][3*f.y_jos+2]=C[0];
    }
    for(i=f.y_sus;i<=f.y_jos;i++)
    {
        I[f.x_sus][3*i]=I[f.x_jos][3*i]=C[2];
        I[f.x_sus][3*i+1]=I[f.x_jos][3*i+1]=C[1];
        I[f.x_sus][3*i+2]=I[f.x_jos][3*i+2]=C[0];
    }

    fseek(fin,54,SEEK_SET);

    unsigned char completare[3]={0};

    for(i=h_img-1;i>=0;i--)
    {
        fwrite(I[i],3*w_img,1,fin);
        if(padding)
            fwrite(completare,padding,1,fin);
    }

    fclose(fin);

    for(i=0;i<h_img;i++)
        free(I[i]);
    free(I);
}

int cmp(void const *a, void const *b)
{
    detectie p=*(detectie *)a,t=*(detectie *)b;

    if(p.corelatie<t.corelatie)return 1;
    if(p.corelatie>t.corelatie)return -1;
    return 0;
}

void sortare(detectie *D,int n)
{
    qsort(D,n,sizeof(detectie),cmp);
}

float suprapunere(fereastra a,fereastra b)
{
    int h_a,w_a,h_b,w_b;

    h_a=a.x_jos-a.x_sus+1;
    w_a=a.y_jos-a.y_sus+1;

    h_b=b.x_jos-b.x_sus+1;
    w_b=b.y_jos-b.y_sus+1;

    int h_intersectie,w_intersectie;

    if(a.x_sus<=b.x_sus)
       {
          if(b.x_sus-a.x_sus+1>h_a)return 0;
          h_intersectie=a.x_jos-b.x_sus+1;

          if(a.y_sus<=b.y_sus)
           {
              if(b.y_sus-a.y_sus+1>w_a)return 0;
              w_intersectie=a.y_jos-b.y_sus+1;
           }
          else
           {
              if(a.y_sus-b.y_sus+1>w_b)return 0;
              w_intersectie=b.y_jos-a.y_sus+1;
           }
       }
    else
       {
          if(a.x_sus-b.x_sus+1>h_b)return 0;
          h_intersectie=b.x_jos-a.x_sus+1;

          if(b.y_sus<=a.y_sus)
           {
              if(a.y_sus-b.y_sus+1>w_b)return 0;
              w_intersectie=b.y_jos-a.y_sus+1;
           }
          else
           {
              if(b.y_sus-a.y_sus+1>w_a)return 0;
              w_intersectie=a.y_jos-b.y_sus+1;
           }
       }

    float arie_a,arie_b,arie_intersectie;

    arie_a=w_a*h_a;
    arie_b=w_b*h_b;
    arie_intersectie=w_intersectie*h_intersectie;

    float raport;

    raport=arie_intersectie/(arie_a+arie_b-arie_intersectie);

    return raport;

}

void eliminare_nonmaxime(detectie **D,int *nr_detectii)
{
      sortare((*D),*nr_detectii);

      int i,j,p;

      for(i=0;i<*nr_detectii-1;i++)
        for(j=i+1;j<*nr_detectii;j++)
         if(suprapunere((*D)[i].punct,(*D)[j].punct)>0.2)
            {
                for(p=j;p<*nr_detectii-1;p++)
                (*D)[p]=(*D)[p+1];

                (*nr_detectii)--;

                *D=(detectie *)realloc(*D,(*nr_detectii)*sizeof(detectie));

                if(*D==NULL)
                {
                    printf("Eroare la realocarea memoriei pentru vectorul de detectii in cadrul stergerii!\n");
                    return;
                }

                j--;
            }
}

int main()
{
    char *nume_imagine,**sablon,*nume_fisier_in;
    unsigned char **culoare;

    nume_fisier_in=(char *)malloc(100*sizeof(char));
    if(nume_fisier_in==NULL)
    {
        printf("Nu am putut aloca memorie pentru tabloul ce retine numele fisierului care contine datele de intrare!\n");
        return 0;
    }

    printf("Numele fisierului din care citesc numele imaginii, numele sabloanelor si culorile sabloanelor este: \n");

    fgets(nume_fisier_in,100,stdin);

    nume_fisier_in[strlen(nume_fisier_in)-1]='\0';

    FILE *f=fopen(nume_fisier_in,"r");

    if(f==NULL)
    {
        printf("Eroare la deschiderea fisierului care contine datele de intrare!\n");
        return 0;
    }

    nume_imagine=(char *)malloc(100*sizeof(char));
    if(nume_imagine==NULL)
    {
        printf("Eroare la alocarea memoriei pentru tabloul care retine calea imaginii!\n");
        return 0;
    }

    fgets(nume_imagine,100,f);

    nume_imagine[strlen(nume_imagine)-1]='\0';

    char nume_imagine_gs[]="imagine_gri.bmp";

    grayscale_image(nume_imagine,nume_imagine_gs);

    sablon=NULL;
    culoare=NULL;
    int nr_sabloane=0,i;

    while(!feof(f))
    {
        sablon=(char **)realloc(sablon,(nr_sabloane+1)*sizeof(char *));

        if(sablon==NULL)
        {
            printf("Eroare la realocarea memoriei pentru tabloul care memoreaza caile sabloanelor!\n");
            return 0;
        }

        culoare=(unsigned char **)realloc(culoare,(nr_sabloane+1)*sizeof(unsigned char *));

        if(culoare==NULL)
        {
            printf("Eroare la realocarea memoriei pentru tabloul care memoreaza culorile sabloanelor!\n");
            return 0;
        }

        sablon[nr_sabloane]=(char *)malloc((L_tmpnam+1)*sizeof(char));

        if(sablon[nr_sabloane]==NULL)
        {
            printf("Eroare la alocarea memoriei pentru sablonul %d \n",nr_sabloane);
            return 0;
        }

        culoare[nr_sabloane]=(unsigned char *)malloc(3*sizeof(unsigned char));

        if(culoare==NULL)
        {
            printf("Eroare la alocarea memoriei pentru culoarea sablonului %d !\n",nr_sabloane);
            return 0;
        }

        char *str;
        str=(char *)malloc(100*sizeof(char));

        if(str==NULL)
        {
            printf("Eroare la alocarea memoriei pentru tabloul in care retin pe rand numele fisierelor!\n");
            return 0;
        }

        tmpnam(sablon[nr_sabloane]);

        fscanf(f,"%s%u%u%u\n",str,&culoare[nr_sabloane][0],&culoare[nr_sabloane][1],&culoare[nr_sabloane][2]);

        grayscale_image(str,sablon[nr_sabloane]);

        nr_sabloane++;
    }

    detectie *D;
    int nr_detectii;

    nr_detectii=0;
    D=NULL;


    for(i=0;i<nr_sabloane;i++)
       template_matching(nume_imagine_gs,sablon[i],0.5,&D,&nr_detectii,i);

    if(D!=NULL)
       eliminare_nonmaxime(&D,&nr_detectii);

    if(D!=NULL)
      for(i=0;i<nr_detectii;i++)
         colorare(nume_imagine,D[i].punct,culoare[D[i].nr_culoare]);

    fclose(f);
    free(D);
    for(i=0;i<nr_sabloane;i++)
    {
        free(sablon[i]);
        free(culoare[i]);
    }
    free(sablon);
    free(culoare);
    free(nume_imagine);
    free(nume_fisier_in);


    return 0;
}
