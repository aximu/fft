#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#define N 2048
#define PI atan(1)* 4           //pi�ı�ʾ

/*���帴������*/
typedef struct {
    float real;
    float img;
}complex;

complex x[N], * Wn;                     /*��������,����ϵ��Wn*/
//int size_x;                           /*size_xΪ�������źŸ���������Ϊ2�ı���*/
void fft(complex* IN_X, int n);         /*������IN_X�����ٸ���Ҷ�任*/
void InitWn(complex* w, int n);          /*����һ������Ϊn��Wnŷ����ʽ����*/
void Reverse(complex* IN_X, int n);      /*������IN_X��ַ*/
complex add(complex, complex);          /*�����ӷ�*/
complex mul(complex, complex);          /*�����˷�*/
complex sub(complex, complex);          /*��������*/
void output();                          /*������ٸ���Ҷ�任�Ľ��*/


int main()
{
    int i;                             /*������*/
    float dataA[N];

    /*˳������������*/
    for (i = 0; i < N; i++)
    {
        x[i].real = 0.01*i;
        x[i].img = 0;
    //    printf("x[%d].real = %f\n", i, x[i].real);
    }
    /*������ļ�д��*/
    //FILE* fp;
    //float fft_in;
    //if ((fp = fopen("sample_input.txt", "r")) == NULL)
    //{
    //    printf("�޷����ļ�");
    //    exit(0);//��ֹ���� 
    //}
    //for (i = 0; i < N; i++)
    //{
    //    fscanf(fp, "%f", &fft_in);
    //    x[i].real = fft_in;
    //    x[i].img = 0;			//���벨�ε���������
    //    printf("x[%d].real = %f\n", i, fft_in);
    //}
    //fclose(fp);

    Wn = (complex*)malloc(sizeof(complex) * N);  //����任���ֵ���ڴ�ռ�
    if (Wn == NULL)
    {
        printf("Cannot allocate memory for Wn!\n");
        return -1;
    }
    
    InitWn(Wn, N);           //���ñ任��
    fft(x, N);               //���ÿ��ٸ���Ҷ�任
    printf("���FFT��Ľ��\n");
    output();                //�����������Ҷ�任�������
 
    free(Wn);

    //�������Ƶ����,dataA[i] = sqrt(dataR[i]*dataR[i]+dataI[i]*dataI[i])
    for (i = 0; i < N; i++)
    {
        dataA[i] = sqrt((x[i].real) *  (x[i].real) +  (x[i].img) *  (x[i].img));
        dataA[i] = dataA[i] / N;
        printf("dataA[%d]: %f\n", i, dataA[i]);
    }

    return 0;
}


/*���ٸ���Ҷ�任*/
void fft(complex* IN_X, int n)
{
    int row = log(n) / log(2);                   //rowΪFFT��N�������Ӧ���������
    int i = 0, j = 0, k = 0, L = 0;
    complex top, bottom, product;
    Reverse(IN_X, n);                           //���ñ�ַ����
    for (i = 0; i < row; i++)                    /*��i���������� */
    {
        L = 1 << i;                         //L����2��i�η�,��ʾÿһ�����������ľ��룬��һ���е��θ���
        for (j = 0; j < n; j += 2 * L)       /*jΪ��ƫ�Ƶ�ַ��ÿL��������һ�飬ÿ����N/2L��*/
        {
            for (k = 0; k < L; k++)          /*kΪһ���е��ε�ƫ�Ƶ�ַ��һ���������� ÿ��group�ڵĵ�������*/
            {
                product = mul(IN_X[j + k + L], Wn[n * k / 2 / L]);
                top = add(IN_X[j + k], product);
                bottom = sub(IN_X[j + k], product);
                IN_X[j + k] = top;
                IN_X[j + k + L] = bottom;

                //IN_X[j + k].real = IN_X[j + k].real/2;
                //IN_X[j + k].img = IN_X[j + k].img / 2;
                //IN_X[j + k + L].real = IN_X[j + k + L].real/ 2;
                //IN_X[j + k + L].img = IN_X[j + k + L].img / 2;

            //    printf("product(%d %d %d) = %f + j*%f \n", i, j, k,  (product.real),  (product.img));
            //   printf("mul1: %f %f\t mul2: %f %f\n",  (IN_X[j + k + L].real),  (IN_X[j + k + L].img), Wn[n * k / 2 / L].real, Wn[n * k / 2 / L].img);
            }
        }
    }
}


/*����һ������ŷ����ʽ��Wn��ֵ*/
void InitWn(complex* w, int n)   //nΪ����ĸ���,wΪȨֵWn
{
    int k;
    for (k = 0; k < n; k++)
    {
        w[k].real = cos(2 * PI / n * k);   //��ŷ����ʽ������ת����
        w[k].img = -1 * sin(2 * PI / n * k);
     //   printf("w[%d].real = %f\tw[%d].img = %f\n", k, w[k].real,  k, w[k].img);
    }
}


/*��ַ���㣬��x(n)��λ����*/
void Reverse(complex* IN_X, int n)
{
    int row = log(n) / log(2);  //rowΪFFT��N�������Ӧ���������
    complex temp;               //��ʱ��������
    unsigned short i = 0, j = 0, k = 0;
    unsigned int m, t, F0, F1;
 
    for (i = 0; i < N; i++)								 //������Ԫ��ִ����䵹��
    {
        /*��ȡ�±�I�ķ���J����ֵ*/
        j = 0;
        for (k = 0; k < (row / 2 + 0.5); k++)					//k��ʾ����
        {
            //*�������*/
            m = 1;									//m�����λΪ1�Ķ�������
 //           t = (unsigned int)pow(2, row - 1);			//t�ǵ�MλΪ1�Ķ�������
            t = 1 << (row - 1);
            m <<= k;								//��m����kλ
            t >>= k;								//��n����kλ
            F0 = i & t;								//I��n��λ����ȡ��ǰ�벿�ֵ�kλ
            F1 = i & m;								//I��m��λ����ȡ��F0��Ӧ�ĺ�벿�ֵĵ�λ
            if (F0) j = j | m;						//J��m��λ��ʹF0��Ӧ��λΪ1
            if (F1) j = j | t;						//J��n��λ��ʹF1��Ӧ��λΪ1 
        }

        if (i < j)
        {
            temp = IN_X[i];
            IN_X[i] = IN_X[j];
            IN_X[j]= temp;
        }
    }
}

/*�������Ҷ�任�Ľ��*/
void output()
{
    int i;
    printf("The result are as follows��\n");
    for (i = 0; i < N; i++)
    {
        printf("%.4f\t", x[i].real);
        if (x[i].img >= 0.0001)
            printf("+%.5fj\n", x[i].img);
        else if (fabs(x[i].img) < 0.0001)
            printf("\n");
        else
            printf("%.5fj\n", x[i].img);
    }
}

//�����ӷ��Ķ���
complex add(complex a, complex b)
{
    complex c;
    c.real = a.real + b.real;
    c.img = a.img + b.img;
    return c;
}

//�����˷��Ķ���
complex mul(complex a, complex b)
{
    complex c;
    c.real = a.real * b.real - a.img * b.img;
    c.img = a.real * b.img + a.img * b.real;
    return c;
}

//���������Ķ���
complex sub(complex a, complex b)
{
    complex c;
    c.real = a.real - b.real;
    c.img = a.img - b.img;
    return c;
}

