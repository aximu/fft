#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#define N 2048
#define PI atan(1)* 4           //pi的表示

/*定义复数类型*/
typedef struct {
    float real;
    float img;
}complex;

complex x[N], * Wn;                     /*输入序列,复数系数Wn*/
//int size_x;                           /*size_x为采样的信号个数，必须为2的倍数*/
void fft(complex* IN_X, int n);         /*对输入IN_X，快速傅里叶变换*/
void InitWn(complex* w, int n);          /*生成一个长度为n的Wn欧拉形式序列*/
void Reverse(complex* IN_X, int n);      /*对输入IN_X地址*/
complex add(complex, complex);          /*复数加法*/
complex mul(complex, complex);          /*复数乘法*/
complex sub(complex, complex);          /*复数减法*/
void output();                          /*输出快速傅里叶变换的结果*/


int main()
{
    int i;                             /*输出结果*/
    float dataA[N];

    /*顺序生成书输入*/
    for (i = 0; i < N; i++)
    {
        x[i].real = 0.01*i;
        x[i].img = 0;
    //    printf("x[%d].real = %f\n", i, x[i].real);
    }
    /*输入从文件写入*/
    //FILE* fp;
    //float fft_in;
    //if ((fp = fopen("sample_input.txt", "r")) == NULL)
    //{
    //    printf("无法打开文件");
    //    exit(0);//终止程序 
    //}
    //for (i = 0; i < N; i++)
    //{
    //    fscanf(fp, "%f", &fft_in);
    //    x[i].real = fft_in;
    //    x[i].img = 0;			//输入波形的虚数部分
    //    printf("x[%d].real = %f\n", i, fft_in);
    //}
    //fclose(fp);

    Wn = (complex*)malloc(sizeof(complex) * N);  //分配变换后的值的内存空间
    if (Wn == NULL)
    {
        printf("Cannot allocate memory for Wn!\n");
        return -1;
    }
    
    InitWn(Wn, N);           //调用变换核
    fft(x, N);               //调用快速傅里叶变换
    printf("输出FFT后的结果\n");
    output();                //调用输出傅里叶变换结果函数
 
    free(Wn);

    //求出幅度频率谱,dataA[i] = sqrt(dataR[i]*dataR[i]+dataI[i]*dataI[i])
    for (i = 0; i < N; i++)
    {
        dataA[i] = sqrt((x[i].real) *  (x[i].real) +  (x[i].img) *  (x[i].img));
        dataA[i] = dataA[i] / N;
        printf("dataA[%d]: %f\n", i, dataA[i]);
    }

    return 0;
}


/*快速傅里叶变换*/
void fft(complex* IN_X, int n)
{
    int row = log(n) / log(2);                   //row为FFT的N个输入对应的最大列数
    int i = 0, j = 0, k = 0, L = 0;
    complex top, bottom, product;
    Reverse(IN_X, n);                           //调用变址函数
    for (i = 0; i < row; i++)                    /*第i级蝶形运算 */
    {
        L = 1 << i;                         //L等于2的i次方,表示每一级蝶形中组间的距离，即一组中蝶形个数
        for (j = 0; j < n; j += 2 * L)       /*j为组偏移地址，每L个蝶形是一组，每级有N/2L组*/
        {
            for (k = 0; k < L; k++)          /*k为一组中蝶形的偏移地址，一个蝶形运算 每个group内的蝶形运算*/
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


/*产生一个周期欧拉形式的Wn的值*/
void InitWn(complex* w, int n)   //n为输入的个数,w为权值Wn
{
    int k;
    for (k = 0; k < n; k++)
    {
        w[k].real = cos(2 * PI / n * k);   //用欧拉公式计算旋转因子
        w[k].img = -1 * sin(2 * PI / n * k);
     //   printf("w[%d].real = %f\tw[%d].img = %f\n", k, w[k].real,  k, w[k].img);
    }
}


/*变址计算，将x(n)码位倒置*/
void Reverse(complex* IN_X, int n)
{
    int row = log(n) / log(2);  //row为FFT的N个输入对应的最大列数
    complex temp;               //临时交换变量
    unsigned short i = 0, j = 0, k = 0;
    unsigned int m, t, F0, F1;
 
    for (i = 0; i < N; i++)								 //对数组元素执行码间倒序
    {
        /*获取下标I的反序J的数值*/
        j = 0;
        for (k = 0; k < (row / 2 + 0.5); k++)					//k表示操作
        {
            //*反序操作*/
            m = 1;									//m是最低位为1的二进制数
 //           t = (unsigned int)pow(2, row - 1);			//t是第M位为1的二进制数
            t = 1 << (row - 1);
            m <<= k;								//对m左移k位
            t >>= k;								//对n右移k位
            F0 = i & t;								//I与n按位与提取出前半部分第k位
            F1 = i & m;								//I与m按位与提取出F0对应的后半部分的低位
            if (F0) j = j | m;						//J与m按位或使F0对应低位为1
            if (F1) j = j | t;						//J与n按位或使F1对应高位为1 
        }

        if (i < j)
        {
            temp = IN_X[i];
            IN_X[i] = IN_X[j];
            IN_X[j]= temp;
        }
    }
}

/*输出傅里叶变换的结果*/
void output()
{
    int i;
    printf("The result are as follows：\n");
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

//复数加法的定义
complex add(complex a, complex b)
{
    complex c;
    c.real = a.real + b.real;
    c.img = a.img + b.img;
    return c;
}

//复数乘法的定义
complex mul(complex a, complex b)
{
    complex c;
    c.real = a.real * b.real - a.img * b.img;
    c.img = a.real * b.img + a.img * b.real;
    return c;
}

//复数减法的定义
complex sub(complex a, complex b)
{
    complex c;
    c.real = a.real - b.real;
    c.img = a.img - b.img;
    return c;
}

