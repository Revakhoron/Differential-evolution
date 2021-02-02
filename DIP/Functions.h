#pragma once
#pragma once
#include <iostream>
#include <vector>

#define _USE_MATH_DEFINES
#include <math.h>

#define dimension 2

float maxVal = 0;
float minVal = 0;


inline float ackley(std::vector<float> params)
{
    float eval_value = 32.768 * 2;
    for (int i = 0; i < dimension; i++)
    {
        params[i] = (params[i] - 0.5) * eval_value;
    }
    float a = 20.0f, b = 0.2f, c = 2 * M_PI, d = dimension;
    float sum1 = 0, sum2 = 0;
    float s1 = 0, s2 = 0, r_value = 0;
    for (int i = 0; i < d; i++)
    {
        sum1 += pow(params[i], 2);
        sum2 += cos(c * params[i]);
    }
    s1 = -a * exp((-b) * sqrt(sum1 / d));
    s2 = -exp(sum2 / d);
    r_value = s1 + s2 + a + (exp(1));
    return r_value;
}


inline float griewank(std::vector<float> params)
{
    float eval_value = 600 * 2;
    float d = dimension;
    for (int i = 0; i < dimension; i++)
    {
        params[i] = (params[i] - 0.5) * eval_value;
    }
    float sum1 = 0, sum2 = 1;
    for (int i = 1; i <= d; i++)
    {
        sum1 += ((pow(params[i - 1], 2)) / 4000);
        sum2 *= cos(params[i - 1] / sqrt(i));
    }

    return(sum1 - sum2 + 1.0f);
}

inline float rosenbrock(std::vector<float> params)
{
    float eval_value = 2.048 * 2;
    float d = dimension;
    for (int i = 0; i < dimension; i++)
    {
        params[i] = (params[i] - 0.5) * eval_value;
    }
    float sum = 0;
    for (int i = 0; i < d - 1; i++)
    {
        sum += (100 * pow((params[i + 1] - pow(params[i], 2)), 2) + (pow((params[i] - 1), 2)));
    }
    return sum;
}

inline float wi(float xi)
{
    return(1.0 + ((xi - 1.0) / 4.0));
}


inline float levy(std::vector<float> params)
{
    float eval_value = 10.0 * 2;
    float d = dimension;
    for (int i = 0; i < dimension; i++)
    {
        params[i] = (params[i] - 0.5) * eval_value;
    }
    float sum = 0;
    for (int i = 0; i < d - 1; i++)

    {
        sum += ((pow(wi(params[i]) - 1, 2) * (1 + 10 * pow(sin(M_PI * wi(params[i]) + 1), 2))) + (pow(wi(wi(params[d - 1]) - 1), 2) * (1 + pow(sin(2 * M_PI * wi(params[d - 1])), 2))));
    }
    return (pow(sin(M_PI * wi(params[0])), 2) + sum);
}

inline float rastring(std::vector<float> params)
{
    float eval_value = 5.12 * 2;
    float d = dimension;
    for (int i = 0; i < dimension; i++)
    {
        params[i] = (params[i] - 0.5) * eval_value;
    }
    float sum = 0;
    for (int i = 0; i < d; i++)
    {
        sum += (pow(params[i], 2) - 10 * cos(2 * M_PI * params[i]));
    }

    return(10 * d + sum);
}


inline float schwefel(std::vector<float> params)
{
    float sum = 0;
    float eval_value = 500.0 * 2;
    float d = dimension;
    for (int i = 0; i < dimension; i++)
    {
        params[i] = (params[i] - 0.5) * eval_value;
    }
    for (int i = 0; i < d; i++)
    {
        sum += params[i] * sin(sqrt(abs(params[i])));
    }

    return (418.9829 * d - sum);
}

inline float sphere(std::vector<float> params)
{
    float sum = 0;
    float eval_value = 5.12 * 2;
    float d = dimension;
    for (int i = 0; i < dimension; i++)
    {
        params[i] = (params[i] - 0.5) * eval_value;
    }

    for (int i = 0; i < d; i++)
    {
        sum += pow(params[i], 2);
    }
    return sum;
}

inline float zakharov(std::vector<float> params)
{
    float eval_value = 10.0;
    float d = dimension;
    for (int i = 0; i < dimension; i++)
    {
        params[i] = (params[i] - 0.5) * eval_value;
        //if (params[i] > 0)
        //{
          //  params[i] = params[i] * 2;
        //}
    }

    float sum1 = 0, sum2 = 0;

    for (int i = 0; i < d; i++)
    {
        sum1 += pow(params[i], 2);
        sum2 += 0.5 * (i + 1.0) * params[i];
    }

    return (sum1 + pow(sum2, 2) + pow(sum2, 4));
}

inline float michalewicz(std::vector<float> params)
{
    float eval_value = M_PI * 2;
    float d = dimension;
    for (int i = 0; i < dimension; i++)
    {
        params[i] = (params[i]) * eval_value;
    }

    float sum = 0;

    for (int i = 0; i < d; i++)
    {
        sum += sin(params[i]) * pow(sin((float(i) + 1.0) * pow(params[i], 2) / M_PI), 20);
    }
    return (-sum);
}