#include <stdio.h>

int main(){
    printf("single\n");
    int x = 100;
    #pragma omp parallel for firstprivate(x) lastprivate(x)
    for (int i = 1;i < 10;i++) {
        printf("test2, %d\n", x);
    }
    printf("single %d\n", x);
}
