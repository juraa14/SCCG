#include <iostream>
#include <chrono>

int main()
{   

    auto start = std::chrono::high_resolution_clock::now();

    std::cout << "Hello World!\n";

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Execution time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
}

