#include <cstdio>
#include <cstring>
#include <cstdlib>
using namespace std;
class Node
{
public:
    int r,c;
    int *arr;
public:
    Node(int a, int b): r(a), c(b)
    {
        arr = (int*) calloc (a*b, sizeof (int));
    }
    ~Node()
    {
        printf("%d:", arr);
        free(arr); 
        puts("Freed!");
    }
    friend Node ttt();
};

Node ttt()
{
    Node ret(2,2);
    for (int i=0;i<2;i++)
        for (int j=0;j<2;j++) scanf("%d",&ret.arr[i*2 + j]);
    return ret;
}
int main()
{
    Node me = ttt();
    // printf("%d\n",me.arr);
    return 0;
}