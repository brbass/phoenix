#include <iostream>
#include <stdexcept>
#include <string>

void Check(bool condition, std::string message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}
