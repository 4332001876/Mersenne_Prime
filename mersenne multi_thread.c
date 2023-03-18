//罗浩铭 PB21030838 所有代码均为原创,每个字符均为手打 版权所有
#include <stdio.h>
#include <windows.h>
#include <time.h>
#include <pthread.h>
#define LENTH 64
#define DIV 429496729 //(2^32)/10
#define REM 6         //(2^32)%10
#define FILTER 0xffffffffULL

int multi_thread_partition[17] = {3, 341, 461, 569, 619, 661, 709, 751, 787, 821, 839, 863, 887, 929, 953, 977, 1001}; //多线程的任务分配

int is_prime[2010] = {0};          // 若is_prime[i]的值为1代表i是质数,值为0代表i是合数
int is_mersenne_prime[1001] = {0}; // 若is_prime[p]的值为1代表2^p-1是质数,值为0代表2^p-1是合数

void *thread_fun(void *p); //多线程运行的任务函数

int miller_rabin(int base, int p); //返回1为质数,0为非质数
// lucas比miller_rabin快10%左右,而且是充要条件,没有出错概率,所以用lucas,miller_rabin只是个摆设
int lucas(int p); //卢卡斯莱默检验法,返回1为质数,0为非质数

void mersenne_mod(unsigned long long *a, int p);      // (a^(2^p-1)) % (2^p-1),会破坏a,这是整个程序最重要的优化!
void mersenne_mod_unit(unsigned long long *a, int p); // 梅森素数模子操作,(a^(2^p-1)) % (2^p-1),会破坏a
//原理讲解:考虑十进制，想象一下12345678 % 999，同余之后变成45678+12300=57978，再同余变成978+57=1035,再同余变成1+35=36
//再考虑二进制的情况，梅森数是2^p-1，同样地每次截取左边一部分向右移p个二进制位加到右边去,这就是这个算法的原理

void small_prime(int *is_prime);                             //返回0-2009的质数信息(是否为质数)
void clear(unsigned long long *a);                           //将a的每一位数变成0
int mod10(unsigned long long *a);                            //求a%10
int convert10(unsigned long long *a, char *str);             //将a转换为10进制形式的字符串,以便显示,对a是破坏性的,会将a置零
void generate_mersenne(unsigned long long *mersenne, int p); //清空mersenne数组,并在数组里生成梅森数2^p-1
int mod_smallint(unsigned long long *a, int divisor);        //要满足divisor<2^30,不破坏a
int cmp(unsigned long long *a, unsigned long long *b);       //比较大小,a>b返回1, a==b返回0, a<b返回-1

void arraycopy(unsigned long long *dest, unsigned long long *b); // copy到dest里
void display(unsigned long long *a);                             //调试用,输出数组每一位
void display_str(char *str);                                     //打印十进制值的字符串
void show_mersenne(int p);                                       //打印梅森素数的十进制值

int add(unsigned long long *a, unsigned long long *b);                           // a+b
int add_c(unsigned long long *a, unsigned long long *b, unsigned long long *c);  // a+b,结果装c里
int sub(unsigned long long *a, unsigned long long *b);                           // a-b,结果装a里,请保证a>b
int sub_c(unsigned long long *a, unsigned long long *b, unsigned long long *c);  // a-b,结果装c里,请保证a>b
int sub_int(unsigned long long *a, unsigned long long b);                        // a-b,b为小整数,结果装a里,请保证a>b
void multi(unsigned long long *a, unsigned long long *b, unsigned long long *c); // a*b,结果装c里
void multi_int(unsigned long long *a, int num);                                  // a*num,结果装a里

int main(void)
{

    time_t start_t, finish_t;
    int count = 0;
    unsigned long long mersenne[LENTH] = {0};

    char output[10 * LENTH + 1] = {0}; //最后一位封尾

    start_t = clock();

    small_prime(is_prime);
    int i;
    is_mersenne_prime[2] = 1;
    printf("2^2-1    3\n"); // 2为偶数,特判

    pthread_t thread[16];
    int series_num[16] = {0}; //拿来协调多线程的数组
    for (i = 0; i < 16; i++)
    {
        series_num[i] = i;
    }
    for (i = 0; i < 16; i++)
    {
        pthread_create(&thread[i], NULL, thread_fun, &series_num[i]);
    }
    for (i = 0; i < 16; i++)
    {
        pthread_join(thread[i], NULL);
    }
    for (i = 3; i < 1000; i++)
    {
        if (is_mersenne_prime[i])
        {
            generate_mersenne(mersenne, i);
            printf("\n2^%d-1    ", i);

            convert10(mersenne, output);
            display_str(output);
        }
    }
    finish_t = clock();
    printf("\nIt took %dms to complete.", finish_t - start_t); //输出运行时间
    system("pause");
    return 0;
}

void *thread_fun(void *p) //多线程运行的任务函数
{
    int num = *((int *)p);
    int i;
    unsigned long long mersenne[LENTH] = {0};
    char output[10 * LENTH + 1] = {0}; //最后一位封尾
    for (i = multi_thread_partition[num]; i < multi_thread_partition[num + 1]; i += 2)
    {
        generate_mersenne(mersenne, i);
        //一堆筛
        // i必须为素数
        if (!is_prime[i])
            continue;
        //模4余3筛,若p大于3且模4余3,则2p+1必须为合数,否则2^p-1非素数
        if (i > 3 && i % 4 == 3 && is_prime[2 * i + 1])
            continue;
        //整除筛,2^p-1所有素因子可表示为2kp+1的形式，可用此法简单地筛去一些数,这个筛可以让程序快大概10-20%
        if (i > 10)
        {
            if (mod_smallint(mersenne, 2 * i + 1) == 0) //不破坏mersenne
                continue;
            if (mod_smallint(mersenne, 6 * i + 1) == 0) //不破坏mersenne
                continue;
            if (mod_smallint(mersenne, 8 * i + 1) == 0) //不破坏mersenne
                continue;
            if (mod_smallint(mersenne, 10 * i + 1) == 0) //不破坏mersenne
                continue;
            if (mod_smallint(mersenne, 24 * i + 1) == 0) //不破坏mersenne
                continue;
        }
        //筛部分结束,开始正常检验
        // lucas比miller_rabin快10%左右,而且是充要条件,所以用lucas
        if (lucas(i))
        {
            is_mersenne_prime[i] = 1;
        }
    }
}
int lucas(int p) //卢卡斯莱默检验法,返回1为质数,0为非质数
{
    unsigned long long one[LENTH] = {0};
    one[LENTH - 1] = 1;
    unsigned long long temp;
    unsigned long long *a = (unsigned long long *)malloc(LENTH * sizeof(unsigned long long));
    unsigned long long *b = (unsigned long long *)malloc(LENTH * sizeof(unsigned long long));
    unsigned long long *c = (unsigned long long *)malloc(LENTH * sizeof(unsigned long long));
    unsigned long long *temp_p;
    clear(a);
    clear(b);
    clear(c);

    a[LENTH - 1] = 4;
    b[LENTH - 1] = 4;
    int i = 0;
    for (i = 0; i < p - 2; i++)
    {
        multi(a, b, c);
        if (c[LENTH - 1] >= 2) //减2操作
        {
            c[LENTH - 1] -= 2;
        }
        else
        {
            if (cmp(c, one) <= 0) //若c=1,则将c-2转换为p-1;若c=0,则将c-2转换为p-2
            {
                temp = 2 - c[LENTH - 1];
                generate_mersenne(c, p);
                c[LENTH - 1] -= temp;
            }
            else
            {
                sub_int(c, 2);
            }
        }

        mersenne_mod(c, p);
        arraycopy(a, c);
        temp_p = c;
        c = b;
        b = temp_p;
    }

    int flag = 1;
    for (i = 0; i < LENTH; i++)
    {
        if (a[i])
            flag = 0; //非零数
    }
    //则a,b装的是2^(p-1)-1次方模,c装的是2^p-2次方模
    free(a);
    free(b);
    free(c);
    if (flag)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}
int miller_rabin(int base, int p) //返回1为质数,0为非质数
{
    unsigned long long one[LENTH] = {0};
    one[LENTH - 1] = 1;
    unsigned long long p_minus_one[LENTH] = {0};
    generate_mersenne(p_minus_one, p - 1); //注意是p-1!
    p_minus_one[LENTH - 1]--;
    unsigned long long *a = (unsigned long long *)malloc(LENTH * sizeof(unsigned long long));
    unsigned long long *b = (unsigned long long *)malloc(LENTH * sizeof(unsigned long long));
    unsigned long long *c = (unsigned long long *)malloc(LENTH * sizeof(unsigned long long));
    unsigned long long *temp_p;
    clear(a);
    clear(b);
    clear(c);

    a[LENTH - 1] = base;
    b[LENTH - 1] = base;
    int i = 0;
    for (i = 0; i < p - 2; i++)
    {
        multi(a, b, c);
        multi_int(c, base);

        mersenne_mod(c, p);
        arraycopy(a, c);
        temp_p = c;
        c = b;
        b = temp_p;
    }

    multi(a, b, c);
    mersenne_mod(c, p);
    //则a,b装的是2^(p-1)-1次方模,c装的是2^p-2次方模
    if (cmp(one, c) == 0)
    {
        free(a);
        free(b);
        free(c);
        return 1;
    }
    else
    {
        if (cmp(one, a) == 0 || cmp(p_minus_one, a) == 0)
        {
            free(a);
            free(b);
            free(c);
            return 1;
        }
        free(a);
        free(b);
        free(c);
        return 0;
    }
}
void small_prime(int *is_prime) //返回0-2009的质数信息(是否为质数)
{
    int i, j;
    for (i = 2; i < 2010; i++)
    {
        is_prime[i] = 1;
    }

    for (i = 2; i < 50; i++)
    {
        if (is_prime[i])
        {
            for (j = 2; j * i < 2010; j++)
            {
                is_prime[j * i] = 0;
            }
        }
    }
}

void arraycopy(unsigned long long *dest, unsigned long long *b) // copy到dest里
{
    int i;
    for (i = 0; i < LENTH; i++)
    {
        dest[i] = b[i];
    }
}
void display(unsigned long long *a) //调试用,输出数组每一位
{
    int i = 0, flag = 0;
    while (i < LENTH && a[i] == 0)
        i++;
    printf("%lu ", a[i]);
    i++;
    while (i < LENTH)
    {
        printf("%lu ", a[i]); //只截取前32位
        i++;
    }
}
void display_str(char *str) //打印十进制值的字符串
{
    int i = 0, flag = 0;
    while (i < 10 * LENTH && (str[i] == 0))
        i++;

    // printf("%s", str + i);
    puts(str + i);
}
void clear(unsigned long long *a) //将a的每一位数变成0
{
    int i;
    for (i = 0; i < LENTH; i++)
    {
        a[i] = 0;
    }
}
int mod10(unsigned long long *a) //求a%10
{
    int i;
    int rem = 0; // TODO:LENTH很长时记得更换

    for (i = 0; i < LENTH; i++)
    {

        rem += 6 * (a[i] % 10);
    }
    rem = rem % 10;
    return rem;
}
int convert10(unsigned long long *a, char *str) //将a转换为10进制形式的字符串,以便显示,对a是破坏性的,会将a置零
{
    int i;
    unsigned long long rem = 0;
    int pos = 0;
    int str_pos = 10 * LENTH - 1;
    while (1)
    {
        while (a[pos] == 0 && pos < LENTH)
        {
            pos++;
        }
        if (pos >= LENTH)
            break;
        for (i = pos; i < LENTH - 1; i++)
        {
            rem = a[i] % 10;

            a[i] = a[i] / 10;
            a[i + 1] += (rem << 32);
        }
        rem = a[LENTH - 1] % 10;
        a[LENTH - 1] = a[LENTH - 1] / 10;

        str[str_pos] = '0' + rem;
        str_pos--;
    }
}

void mersenne_mod(unsigned long long *a, int p) // (a^(2^p-1)) % (2^p-1),会破坏a,这是整个程序最重要的优化!
{
    //原理讲解:考虑十进制，想象一下12345678 % 999，同余之后变成45678+12300=57978，再同余变成978+57=1035,再同余变成1+35=36
    //再考虑二进制的情况，梅森数是2^p-1，同样地每次截取左边一部分加到右边去,这就是这个算法的原理
    unsigned long long mersenne[LENTH] = {0};
    generate_mersenne(mersenne, p);
    while (cmp(a, mersenne) == 1)
    {
        mersenne_mod_unit(a, p);
    }
    if ((cmp(a, mersenne) == 0))
        clear(a);
}
void mersenne_mod_unit(unsigned long long *a, int p) // 梅森素数模子操作,(a^(2^p-1)) % (2^p-1),会破坏a
{
    //原理讲解:考虑十进制，想象一下12345678 % 999，同余之后变成45678+12300=57978，再同余变成978+57=1035,再同余变成1+35=36
    //再考虑二进制的情况，梅森数是2^p-1，同样地每次截取左边一部分加到右边去,这就是这个算法的原理
    int pos = 0;
    unsigned long long temp = 0;
    while (pos < LENTH && a[pos] == 0)
    {
        pos++;
    }
    int min_pos = pos;

    while (pos < LENTH - 1 - p / 32) //可以正常移位
    {
        while (a[pos])
        {
            temp = a[pos];
            a[pos] = 0;
            a[pos + p / 32] += temp >> (p % 32);                    //移位
            a[pos + p / 32 + 1] += (temp << (64 - (p % 32))) >> 32; //移位
        }
        pos++;
    }
    unsigned long long max = (1 << (p % 32)) - 1;
    while (a[pos] > max)
    {
        temp = a[pos];
        a[pos] = (a[pos] << (64 - p % 32)) >> (64 - p % 32);
        a[LENTH - 1] += temp >> (p % 32);
    }

    for (pos = LENTH - 1; pos >= min_pos; pos--)
    {
        a[pos - 1] += a[pos] >> 32;
        a[pos] &= FILTER;
    }
}

int mod_smallint(unsigned long long *a, int divisor) //要满足divisor<2^30,不破坏a
{
    //防止计算过程中溢出
    const unsigned long long num = (1ULL << 32) % divisor;
    unsigned long long result = 0;
    unsigned long long pow = 0;
    unsigned long long pow_num = 1;
    while (pow < LENTH)
    {
        result = (result + a[LENTH - 1 - pow] * pow_num) % divisor;
        pow_num = (pow_num * num) % divisor;
        pow++;
    }
    return (int)result;
}
int cmp(unsigned long long *a, unsigned long long *b) //比较大小,a>b返回1, a==b返回0, a<b返回-1
{
    int val = 0;
    int i = 0;
    while (i < LENTH && a[i] == b[i])
    {
        i++;
    }
    if (i >= LENTH)
    {
        return 0;
    }
    if (a[i] > b[i])
    {
        return 1;
    }
    if (a[i] < b[i])
    {
        return -1;
    }
}
void generate_mersenne(unsigned long long *mersenne, int p) //清空mersenne数组,并在数组里生成梅森数2^p-1
{
    clear(mersenne);
    int pos;
    for (pos = LENTH - 1; pos > LENTH - 1 - p / 32; pos--)
    {
        mersenne[pos] = FILTER;
    }
    mersenne[pos] = (1 << (p % 32)) - 1;
}
void show_mersenne(int p) //打印梅森素数的十进制值
{
    unsigned long long mersenne[LENTH];
    generate_mersenne(mersenne, p);
    char output[10 * LENTH + 1] = {0}; //最后一位封尾
    convert10(mersenne, output);
    display_str(output);
}
int add(unsigned long long *a, unsigned long long *b) // a+b
{
    int i;
    unsigned long long temp = 0;
    unsigned long long filter = (1ULL << 32) - 1;
    for (i = LENTH - 1; i >= 0; i--)
    {
        a[i] = a[i] + b[i] + temp;
        temp = (a[i] >> 32);
        a[i] &= filter;
    }
}
int add_c(unsigned long long *a, unsigned long long *b, unsigned long long *c) // a+b,结果装c里
{
    int i;
    unsigned long long temp = 0;
    unsigned long long filter = (1ULL << 32) - 1;
    for (i = LENTH - 1; i >= 0; i--)
    {
        c[i] = a[i] + b[i] + temp;
        temp = (c[i] >> 32);
        c[i] &= filter;
    }
}
int sub(unsigned long long *a, unsigned long long *b) // a-b,结果装a里,请保证a>b
{
    int i;
    unsigned long long temp = 0;
    unsigned long long carry = (1ULL << 32); //进位
    for (i = LENTH - 1; i >= 0; i--)
    {
        if (a[i] >= b[i] + temp)
        {
            a[i] = a[i] - b[i] - temp;
            temp = 0;
        }
        else
        {
            a[i] |= carry;
            a[i] = a[i] - b[i] - temp;
            temp = 1;
        }
    }
}
int sub_c(unsigned long long *a, unsigned long long *b, unsigned long long *c) // a-b,结果装c里,请保证a>b
{
    int i;
    unsigned long long temp = 0;
    unsigned long long carry = (1ULL << 32); //进位
    for (i = LENTH - 1; i >= 0; i--)
    {
        if (a[i] >= b[i] + temp)
        {
            c[i] = a[i] - b[i] - temp;
            temp = 0;
        }
        else
        {

            c[i] = (a[i] | carry) - b[i] - temp;
            temp = 1;
        }
    }
}
int sub_int(unsigned long long *a, unsigned long long b) // a-b,b为小整数,结果装a里,请保证a>b
{
    int i;
    unsigned long long temp = 0;
    unsigned long long carry = (1ULL << 32); //进位
    for (i = LENTH - 1; i >= 0; i--)
    {
        if (a[i] >= b + temp)
        {
            a[i] = a[i] - b - temp;
            temp = 0;

            break;
        }
        else
        {
            a[i] |= carry;
            a[i] = a[i] - b - temp;
            b = 0;
            temp = 1;
        }
    }
}
void multi(unsigned long long *a, unsigned long long *b, unsigned long long *c) // a*b,结果装c里
{
    clear(c); //清空c,非常重要
    int pos = 0, pos1 = 0, pos2 = 0;

    while (pos1 < LENTH && a[pos1] == 0)
    {
        pos1++;
    }
    int min_pos1 = pos1;
    while (pos2 < LENTH && b[pos2] == 0)
    {
        pos2++;
    }
    int min_pos2 = pos2;
    for (pos1 = LENTH - 1; pos1 >= min_pos1; pos1--)
    {
        for (pos2 = LENTH - 1; pos2 >= min_pos2; pos2--)
        {
            c[-LENTH + 1 + pos1 + pos2] += a[pos1] * b[pos2];
        }
        for (pos = LENTH - 1; pos > min_pos1 + min_pos2 - LENTH; pos--) //请保证LENTH长度预留两倍数字的长度
        {
            c[pos - 1] += c[pos] >> 32;
            c[pos] &= FILTER;
        } //正好不会溢出
    }
}
void multi_int(unsigned long long *a, int num) // a*num,结果装a里
{
    int pos;
    for (pos = 0; pos < LENTH; pos++)
    {
        a[pos] *= num;
    }
    for (pos = LENTH - 1; pos > 0; pos--)
    {
        a[pos - 1] += a[pos] >> 32;
        a[pos] &= FILTER;
    }
}