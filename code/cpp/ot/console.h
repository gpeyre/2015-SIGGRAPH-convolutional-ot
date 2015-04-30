#ifndef CONSOLE_COLOR_H
#define CONSOLE_COLOR_H

#include <iostream>

#if defined(WIN32)
#include <windows.h>
#endif

#define COLOR_WHITE  0
#define COLOR_RED    1
#define COLOR_BLUE   2
#define COLOR_GREEN  3
#define COLOR_YELLOW 4

inline 
std::ostream& blue(std::ostream &s)
{
#if defined(WIN32)
    HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE); 
    SetConsoleTextAttribute(hStdout, FOREGROUND_BLUE|FOREGROUND_GREEN|FOREGROUND_INTENSITY);
#else
    s << "\033[0;34m";
#endif
    return s;
}

inline 
std::ostream& red(std::ostream &s)
{
#if defined(WIN32)
    HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE); 
    SetConsoleTextAttribute(hStdout, FOREGROUND_RED|FOREGROUND_INTENSITY);
#else
    s << "\033[0;31m";
#endif
    return s;
}

inline 
std::ostream& green(std::ostream &s)
{
#if defined(WIN32)
    HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE); 
    SetConsoleTextAttribute(hStdout, FOREGROUND_GREEN|FOREGROUND_INTENSITY);
#else
    s << "\033[0;32m";
#endif
    return s;
}

inline 
std::ostream& yellow(std::ostream &s)
{
#if defined(WIN32)
    HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE); 
    SetConsoleTextAttribute(hStdout, FOREGROUND_GREEN|FOREGROUND_RED|FOREGROUND_INTENSITY);
#else
    s << "\033[0;33m";
#endif
    return s;
}

inline 
std::ostream& white(std::ostream &s)
{
#if defined(WIN32)
    HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE); 
    SetConsoleTextAttribute(hStdout, FOREGROUND_RED|FOREGROUND_GREEN|FOREGROUND_BLUE);
#else
    s << "\033[0;37m";
#endif
    return s;
}

#endif // CONSOLE_COLOR_H
