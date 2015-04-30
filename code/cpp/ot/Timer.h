#ifndef TIMER_H
#define TIMER_H

#include <ctime>
#include <vector>
#include <string>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "console.h"

class Timer
{
public:
    static
    void start(std::vector<double>& timer,
               const int color,
               const std::string msg)
    {
        set_color(color);
        std::cout << msg << white << " ... " << std::endl;
        #ifdef _OPENMP
            timer.push_back(omp_get_wtime());
        #else
            timer.push_back(clock());
        #endif
    }
    
    static
    void start(std::vector<double>& timer,
               const int color,
               const std::string msg,
               const double val)
    {
        set_color(color);
        std::cout << msg << white;
        std::cout << " (" << val << ") ... " << std::endl;
        #ifdef _OPENMP
            timer.push_back(omp_get_wtime());
        #else
            timer.push_back(clock());
        #endif
    }

    static
    void stop(std::vector<double>& timer,
              const int color,
              const std::string msg = std::string())
    {
        double duration = time_duration(timer.back());
        timer.pop_back();
        set_color(color);
        std::cout << "done" << white
        << " (" << duration << " s) " << msg << std::endl;
    }
    
    static
    void set_color(const int color)
    {
        switch (color)
        {
            case COLOR_WHITE:
                std::cout << white;
                break;
            case COLOR_RED:
                std::cout << red;
                break;
            case COLOR_BLUE:
                std::cout << blue;
                break;
            case COLOR_GREEN:
                std::cout << green;
                break;
            case COLOR_YELLOW:
                std::cout << yellow;
                break;
        }
    }
    
    static
    double time_duration(const double init)
    {
        #ifdef _OPENMP
            return (omp_get_wtime() - init);
        #else
            return (clock() - init) / CLOCKS_PER_SEC;
        #endif
    }
};

#endif
