#include <iostream>

// #define DIO_STATS

#include "diophantineSystem.hh"

int main() {
    DiophantineSystem system = DiophantineSystem(6, 6);
    system.insertRow(1, 14, 14);  // 14 = actual sum of row
    system.insertRow(2, 15, 15);  // 15 = actual sum of row
    system.insertRow(2, 17, 17);  // 17 = actual sum of row
    system.insertRow(2, 18, 18);  // 18 = actual sum of row
    system.insertRow(1, 34, 34);  // 34 = actual sum of row
    system.insertRow(2, 15, 15);  // 15 = actual sum of row
    system.insertColumn(26);
    system.insertColumn(28);
    system.insertColumn(32);
    system.insertColumn(25);
    system.insertColumn(41);
    system.insertColumn(26);

    // std::cout << "Solve: " << system.solve() << std::endl;

    int solution_count = 8;
    for(int s = 0; s < solution_count; s++) {
        system.solve();
        // while(system.solve()) {
        std::cout << "\nSolution:" << std::endl;
        for(int row = 0; row < 6; row++) {
            for(int col = 0; col < 6; col++) {
                std::cout << system.solution(row, col) << "  ";
            }
            std::cout << std::endl;
        }
        // system.dumpInfo();
    }

    std::cout << "\nDone!" << std::endl;
    return 0;
}
