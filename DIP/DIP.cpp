// DIP.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/core.hpp>
#include <random>
#include "Functions.h"//import knihoven


//konstanty SF = mutace, CR = krizeni
#define SF 0.5 
#define CR 0.1


//vyber funcke
#define function ackley

//definice barev
#define black cv::Scalar(0, 0, 0) 
#define gray cv::Scalar(150,150,150)
#define red cv::Scalar(0,0,200)
#define blue cv::Scalar(200,100,50)

//velikost okna
int ws_x = 1000;
int ws_y = 1000;


//vykresleni kruhu v danem bode s danou barvou
void MyFilledCircle(cv::Mat& img, cv::Point center, cv::Scalar color)
{
    circle(img,
        center,
        3,
        color,
        cv::FILLED,
        cv::LINE_8);
}

// funkce ktera navraci nejvyssi moznou hodnotu funkce, pro vykresleni pozadi 
template<typename F>//predavam funkci jako argument pomoci template
void testRanges(cv::Mat& img, F func, float& maxVal)
{
    std::vector<float> vec(2);
    vec.emplace_back(0);
    vec.emplace_back(0);
    float curVal;
    for (float y = 0; y < ws_y; y++)
    {
        for (float x = 0; x < ws_x; x++)//iteruji pres jednotlive pixely a zjistuji v jake hodnote ma funkce nejvetsi hodnotu
        {
            vec[0] = x / ws_x;
            vec[1] = y / ws_y;
            curVal = func(vec);
            if (curVal > maxVal)
            {
                maxVal = curVal;
            }
        }
    }
}


template<typename F>
static void colorBackground(cv::Mat& img, float& maxVal, F func)//funkce ktera obarvi pozadi pro urychleni nebarvim po pixelu ale po obdelniku 5x5
{
    std::vector<float> vec(2);
    vec.emplace_back(0);
    vec.emplace_back(0);
    float value;
    for (float y = 0; y < ws_y; y += 5)
    {
        for (float x = 0; x < ws_x; x += 5)
        {
            vec[0] = x / ws_x;
            vec[1] = y / ws_y;
            value = func(vec) / maxVal;
            cv::rectangle(img, cv::Point(x, y), cv::Point(x + 5, y + 5), cv::Scalar(255.0 * value, 255.0 * value, 255.0 * value), -1);
        }
    }
}


struct individual//struktura ktera popisuje jedince jeho pozici a fitness
{
    std::vector<float> position;
    float fitness = 0;
};


inline int rngGen(const int& a, const int& b)//funkce pro generaci nahodnych cisel v rozmezi a-b pro cele cisla
{
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<int> dist6(a, b);

    return dist6(rng);
}


inline float rngGen(const float& a, const float& b)//funkce pro generaci nahodnych cisel v rozmezi a-b pro float cisla
{
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<float> dist6(a, b);

    return (dist6(rng));
}

void generatePopulation(std::vector<individual>& population, const int& NP)//funcke pro generovani populace (nahodna pozice)
{
    population.reserve(NP);//reservuji misto a tim omezim pocet kopirovani
    for (int i = 0; i < NP; i++)//naplnim populaci prazdnymi jedinci
    {
        population.emplace_back(individual());
        std::cout << "emplace_back individual" << std::endl;
    }

    for (int i = 0; i < NP; i++)
    {
        population[i].position.reserve(dimension);
        population[i].position.emplace_back(rngGen(0.0f, 1.0f));
        population[i].position.emplace_back(rngGen(0.0f, 1.0f));//nahodna pozice
    }
}

//funkce pro vykresleni populace
void drawPopulation(std::vector<individual>& population, const int& NP, cv::Mat& image, int& bestIndex)
{
    for (int i = 0; i < NP; i++)
    {
        if (i == bestIndex)
        {//nejlepsi jedinec oznacen jinou barvou
            MyFilledCircle(image, cv::Point(population[i].position[0] * float(ws_x), population[i].position[1] * float(ws_x)), red);
            continue;
        }
        //projdu vsechny jedince a na jejich pozici vykreslim kruh
        MyFilledCircle(image, cv::Point(population[i].position[0] * float(ws_x), population[i].position[1] * float(ws_x)), blue);
    }
    // cv::waitKey(1);
}

//funkce projde vsechny jedince a vypocita se jejich puvodni fitness
template<typename F>
void calculateFitness(std::vector<individual>& pop, const int NP, F func)
{
    for (int i = 0; i < NP; i++)
    {
        pop[i].fitness = func(pop[i].position);
        //pop[i].fitness = ackley(pop[i].position);
    }
}


int getBestIndividualIndex(std::vector<individual>& pop, const int NP, float& bestFitness, int& bestIndex)//funkce navrati nejlepsiho jedince dane populace(respektive jeho index ve vektoru)
{
    for (int i = 0; i < NP; i++)
    {
        if (bestFitness > pop[i].fitness)
        {
            bestFitness = pop[i].fitness;
            bestIndex = i;
        }

    }
    return bestIndex;
}


void checkBounds(individual& mutant)//funkce pro kotrnolu mezi
{
    if (mutant.position[0] > 1.0)
    {
        mutant.position[0] = 1.0;
    }
    if (mutant.position[0] < 0.0)
    {
        mutant.position[0] = 0.0;
    }
    if (mutant.position[1] > 1.0)
    {
        mutant.position[1] = 1.0;
    }
    if (mutant.position[1] < 0.0)
    {
        mutant.position[1] = 0.0;
    }
}

void threeRandom(std::vector<individual>& tr, std::vector<individual>& population, const int& NP, int index)//tri nahodni jedinci
{
    std::vector<int> randoms(3);//pripravim si vektor pro uchovani indexu jedincu
beginAgain:
    int rng;
    int i = 0;
    while (i < 3)
    {
        rng = rngGen(0, NP - 1);
        if (rng != index)
        {
            randoms[i] = rng;
            i++;
        }
    }

    auto it = std::unique(randoms.begin(), randoms.end());//kontroluji unikatnost jedincu ve vektoru
    bool wasUnique = (it == randoms.end());//vysledek ulozim do promenne wasunique

    if (!wasUnique)//pokud nebylo unikatni generuj znova v sekci beginAgain
        goto beginAgain;

    for (int i = 0; i < 3; i++)
    {
        tr[i] = population[randoms[i]];
    }

}

individual mutation(std::vector<individual>& tr)//provedeni mutace ze tri nahodne vybranych
{
    struct individual mutant = {};
    mutant.position.reserve(dimension);
    for (int i = 0; i < dimension; i++)
    {
        mutant.position.emplace_back(0.0f);
    }
    mutant.position[0] = tr[0].position[0] + SF * (tr[1].position[0] - tr[2].position[0]);
    mutant.position[1] = tr[0].position[1] + SF * (tr[1].position[1] - tr[2].position[1]);

    checkBounds(mutant);
    return mutant;
}


int main()
{
    cv::Mat image(ws_x, ws_y + 300, CV_8UC3, cv::Scalar(0, 0, 0));//vytvoreni obrazu pro kresleni
    cv::Mat backGround;//matice ktera bude uchovavat pozadi
    const int NP = 50;//pocet jedincu populace
    std::vector<individual> currentPop;
    std::vector<individual> threeRand(3);
    std::vector<individual> currentPopU(NP);
    std::vector<individual> nextPop(NP);//vektory uchovavajici populace
    individual mutant;//jedinec na ktereho bude uplatnena mutace
    testRanges(image, function, maxVal);

    colorBackground(image, maxVal, function);
    backGround = image.clone();//ulozeni pozadi do opencv matice
    generatePopulation(currentPop, NP);//generace populace
    currentPopU = currentPop;
    nextPop = currentPop;
    int iters = 0;
    int bestIndex = -1;
    float bestFitness = 999999.0;
    cv::imshow("Display window", image);
    int maxIters = 1000;
    while (iters < maxIters)
    {
        image = backGround.clone();//zachovani pozadi
        calculateFitness(currentPop, NP, function);//vypocet fitness u vsech jedincu
        bestIndex = getBestIndividualIndex(currentPop, NP, bestFitness, bestIndex);//index nejlepsiho jedince
        drawPopulation(currentPop, NP, image, bestIndex);//vykresleni populace
        std::stringstream ss;//string pro obraz
        ss << "best fitness: " << currentPop[bestIndex].fitness;
        cv::putText(image, //cilovy obraz
            ss.str(), //text
            cv::Point(ws_y + 25, 15), //pozice
            cv::FONT_HERSHEY_DUPLEX,
            0.5,
            CV_RGB(255, 255, 255), //barva textu
            1);
        ss.str("");
        ss << "Iterace: " << iters;
        cv::putText(image,
            ss.str(),
            cv::Point(ws_y + 25, 30),
            cv::FONT_HERSHEY_DUPLEX,
            0.5,
            CV_RGB(255, 255, 255),
            1);
        ss.str("");
        ss << "x,y: " << currentPop[bestIndex].position[0] << " " << currentPop[bestIndex].position[1];
        cv::putText(image,
            ss.str(),
            cv::Point(ws_y + 25, 45),
            cv::FONT_HERSHEY_DUPLEX,
            0.5,
            CV_RGB(255, 255, 255),
            1);
        cv::imshow("Display window", image);//zobrazeni
        cv::waitKey(33);//cekej 33ms

        for (int i = 0; i < NP; i++)//iteruji pres pocet jedincu
        {
            threeRandom(threeRand, currentPop, NP, i);//tri nahodni jecinci
            mutant = mutation(threeRand);//mutant
            if (rngGen(0.0f, 1.0f) < CR)//pokud je nahodne cislo mensi nez crossover rate tak nahrad jedince mutantem
            {
                currentPopU[i] = mutant;
            }
            else
            {
                currentPopU[i] = currentPop[i];
            }

        }
        calculateFitness(currentPopU, NP, function);//fitness pro trial vektor

        //std::cout << "best fitness is:" << currentPop[bestIndex].fitness << std::endl;
        for (int i = 0; i < NP; i++)
        {
            if (currentPop[i].fitness >= currentPopU[i].fitness)//nahradim stavajici jedince lepsimi nebo stejne dobrymi
            {
                nextPop[i] = currentPopU[i];
            }
            else
            {
                nextPop[i] = currentPop[i];
            }
        }
        iters++;//navysim pocet iteraci
        std::cout << iters << std::endl;
        currentPop = nextPop;
    }
    cv::waitKey(0);
    return 0;
}