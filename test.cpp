#include <cstdio>
using namespace std;
class Node
{
private:
    int r=1,c=2;
    int arr[10];
public:
    Node(int a, int b);
    ~Node();
    // int foo(int x)
    // {
    //     int arr[x];
    //     arr[0]=1;
    // }
};

Node::Node(int a, int b): r(a), c(b)
{
}

Node::~Node()
{
}

int main()
{
    Node X(10, 10);
    printf("%d\n", sizeof (Node));
    return 0;
}