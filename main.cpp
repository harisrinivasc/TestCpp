//
//  main.cpp
//  TestCPP
//
//  Created by Hari Chandrasekar on 6/21/18.
//  Copyright © 2018 Hari Chandrasekar. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <algorithm>
#include <map>
#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <queue>
#include <stack>
#include <string>
#include <random>
#include <bitset>
#include <iomanip>
#include <tuple>
#include <climits>
//#include <cmath>
//#include "Sales_item.h"


//using std::string;
using namespace std;

//int reused = 42; // reused has global scope
typedef unsigned int uint;

class Book
{
    string Title;
    uint Price;
    uint NumSold;
public:
    Book() = default;
    Book(const string s, const uint price, const uint num) : Title(s), Price(price), NumSold(num) {};
    Book(const Book &copy)
    {
        Title = copy.Title;
        this->Price = copy.Price;
        this->NumSold = copy.NumSold;
    }
    void Change(const string s = "", const uint price = 99, const uint num = 99);
    friend istream& Read(istream&, Book&);
    //friend ostream& Print(ostream&, const Book&);
    string GetTitle(void) const;
    uint GetPrice(void) const;
    uint GetNumSold(void) const;
};

istream& Read(istream&, Book&);
ostream& Print(ostream&, const Book&);

istream& Read(istream& is, Book& Book1)
{
    cout << "Enter book Title, Price, NumSold " << endl;
    getline(is, Book1.Title);
    is >> Book1.Price >> Book1.NumSold;
    return is;
}

ostream& Print(ostream& os, const Book& Book1)
{
    os << "Title " << Book1.GetTitle() << ", Price " << Book1.GetPrice() << ", NumSold " << Book1.GetNumSold() << endl;
    return os;
}

void Book::Change(const string s, const uint price, const uint num)
{
    this->Title = s;
    this->NumSold = num;
    this->Price = price;
}

string Book::GetTitle(void) const
{
    return Title;
}

uint Book::GetPrice(void) const
{
    return Price;
}

uint Book::GetNumSold(void) const { return NumSold; }

int f(void);
int g(void);
int h(void);
int j(void);
bool BinSearch(const int, vector<int> *pvec);
int BinSearchRecursive(const int S, const vector<int> &A, const vector<int>::size_type lo, const vector<int>::size_type hi);
void WriteToQueue(queue<string> &WordQ, const string &Word);
long long int Factorial(const int N);
template <typename T> bool Popstack(stack<T> &S, T &Val);
void PrintSetlist(const list<unordered_set<int> > &setlist);
//bool operator< (const string &s1, const string &s2);
void PrintSortVec(const vector<int> &ivec);
void MyMerge(vector<int> &ivec, vector<int>::size_type lo, vector<int>::size_type hi);
void MyMergeSort(vector<int> &ivec, vector<int>::size_type lo, vector<int>::size_type hi);
void MyQuickSort(vector<int> &ivec, vector<int>::size_type lo, vector<int>::size_type hi);
vector<int>::size_type MyQuickPartition(vector<int> &ivec, vector<int>::size_type lo, vector<int>::size_type hi);
void SwapIvec(vector<int> &ivec, vector<int>::size_type lo, vector<int>::size_type hi);
void ShuffleVec(vector<int> &ivec);
void HeapSwimUp(vector<int> &heaptree, vector<int>::size_type k);
void HeapSwimDown(vector<int> &heaptree, vector<int>::size_type k, vector<int>::size_type vecsize);
void InsertHeapTree(vector<int> &heaptree, int val);
int GetHeapTree(vector<int> &heaptree);

#define MAX_COL 5
#define MAX_ROW 3
//#define RC_TO_INT(R, C) (( ( R < MAX_ROW ? R : MAX_ROW - 1 ) * MAX_COL ) + ( C < MAX_COL ? C : MAX_COL - 1 ))
#define RC_TO_INT(R, C) ( ( ( R % MAX_ROW ) * MAX_COL ) + ( C % MAX_COL ) )

class classA
{
private:
    int temp;
public:
    classA()
    {
        temp = 0;
        cout << "Default constructor" << endl;
    }
    classA(const int temp_input) : temp(temp_input)
    {
        cout << "Constrcutor from int: " << temp_input << endl;
    }
    classA(const classA &class_input) : temp(class_input.temp)
    {
        cout << "Copy constructor: " << class_input.temp << endl;
    }
    void operator= (const classA& c1)
    {
        temp = c1.temp;
        cout << "Assignement: " << c1.temp << endl;
    }
    friend bool operator< (const classA& c1, const classA& c2);
    friend classA operator+ (const classA& c1, const classA& c2);
    friend bool operator== (const classA& c1, const classA& c2);
    friend bool operator!= (const classA& c1, const classA& c2);
    classA& operator+= (const classA& c1);
    virtual ~classA()
    {
        cout << "Destroying object: " << temp << endl;
    }
    int getTemp()
    {
        return temp;
    }
    void setTemp(int value)
    {
        temp = value;
    }
};

classA& classA::operator+= (const classA& c1)
{
    temp += c1.temp;
    return *this;
}
classA operator+ (const classA& c1, const classA& c2)
{
    classA sum;
    sum.temp = c1.temp + c2.temp;
    return sum;
}

bool operator< (const classA& c1, const classA& c2)
{
    return (c1.temp < c2.temp);
}

bool operator== (const classA& c1, const classA& c2)
{
    return (c1.temp == c2.temp);
}

bool operator!= (const classA& c1, const classA& c2)
{
    return !(c1 == c2);
}

/*
 //MxN connected component matrix
 #define GET_NUM(R,C)    (((R)*N)+(C))
 const int M = 7, N = 6;
 vector<pair<int,int>> MatrixNgbrs {{-1,0},{0,-1},{1,0},{0,1}};
 vector<vector<int> > ComponentRoot(M,vector<int>(N,-1));
 void SetRoot(const int R, const int C, const int RootVal);
 const vector<vector<int>> ComponentInputs ={{0,1,1,0,0,0},
 {0,0,0,1,1,1},
 {1,1,1,0,0,0},
 {0,0,0,1,0,0},
 {0,1,0,0,0,1},
 {0,1,0,1,0,0},
 {0,1,1,1,0,0}};
 */

int NumWaysStair(const int Total, unordered_map<int,int> &ResultMap, const vector<int> &Steps);
int NumWaysNsum(const int Sum, const vector<int>::size_type VecIdx, unordered_map<string,int> &ResultMap, const vector<int> &InputVec);

class BstType;
class BstType
{
public:
    int key;
    char charkey;
    string value;
    int size;
    shared_ptr<BstType> left;
    shared_ptr<BstType> right;
    BstType() : key(0), value(""), charkey(' '), size(1), left(nullptr), right(nullptr) {}
    ~BstType()
    {
        cout << "Destruct: Key: " << key << " Value: " << value << " Size: " << size << endl;
    }
};

class GraphType
{
private:
    int NumVertex;
    int NumConnections;
    vector<unordered_set<int> > Connections;
public:
    GraphType(const int N)
    {
        NumVertex = N;
        NumConnections = 0;
        for(int i = 0; i < N; ++i)
            Connections.emplace_back(unordered_set<int> ());
    }
    unordered_set<int> &GetConnections(const int Vertex)
    {
        if (Vertex < NumVertex)
            return Connections[Vertex];
        else
            return Connections[0]; //temp. ToDo - should return empty set
    }
    int GetDegree(const int Vertex)
    {
        if (Vertex < NumVertex)
            return Connections[Vertex].size();
        else
            return -9999;
    }
    void AddConnection(const int v1, const int v2)
    {
        if ((v1 < NumVertex) && (v2 < NumVertex))
        {
            //Add connection to both vertices
            auto ret1 = Connections[v1].insert(v2);
            auto ret2 = Connections[v2].insert(v1);
            //If either of the elements was not inserted (connection already present), do not
            //increment connection (edge) count
            if ((true == ret1.second) || (true == ret2.second))
                ++NumConnections;
        }
    }
    void PrintConnections(void)
    {
        cout << "Print connections:" << endl;
        int cnt = 0;
        for (const auto &i : Connections)
        {
            cout << "Vertex: " << cnt++ << " connected to ";
            for (const auto &j : i)
                cout << j << " ";
            cout << endl;
        }
    }
};

class Dfs : public GraphType
{
private:
    vector<bool> Visited;
    vector<int> ParentVertex;
    const int Source;
    int count;
public:
    Dfs(const int N, const int SourceNode) : GraphType(N), Source(SourceNode)
    {
        for (int i = 0; i < N; ++i)
        {
            Visited.push_back(false);
            ParentVertex.push_back(i);
        }
    }
    void PerformDfs()
    {
        cout << "DFS trace from source: " << Source << endl;
        PerformDfs(Source);
        cout << endl;
    }
    void PerformDfs(const int v)
    {
        cout << v << " ";
        Visited[v] = true;
        ++count;
        for(const auto &i : GetConnections(v))
            if (!Visited[i])
            {
                ParentVertex[i] = v;
                PerformDfs(i);
            }
    }
};

class Bfs : public GraphType
{
private:
    vector<bool> Visited;
    vector<int> ParentVertex;
    const int Source;
    queue<int> VertexToVisit;
    int count;
public:
    Bfs(const int N, const int SourceNode) : GraphType(N), Source(SourceNode)
    {
        for (int i = 0; i < N; ++i)
        {
            Visited.push_back(false);
            ParentVertex.push_back(i);
        }
        VertexToVisit.push(Source);
    }
    void PerformBfs(void)
    {
        cout << "BFS trace from source: " << Source << endl;
        Visited[Source] = true;
        cout << Source << " ";
        while (!VertexToVisit.empty())
        {
            int v = VertexToVisit.front();
            VertexToVisit.pop();
            for (const auto &i : GetConnections(v))
                if (!Visited[i])
                {
                    Visited[i] = true;
                    ParentVertex[i] = v;
                    ++count;
                    VertexToVisit.push(i);
                    cout << i << " ";
                }
        }
        cout << endl;
    }
};

Dfs MyDfs(11, 0);
Bfs MyBfs(11, 0);

shared_ptr<BstType> BstHead = nullptr;
void InsertKey(const int key, const string value);
shared_ptr<BstType> InsertKey(shared_ptr<BstType> Node, const int key, const string value);
void DeleteKey(const int key);
shared_ptr<BstType> DeleteKey(shared_ptr<BstType> Node, const int key);
string GetValue(const int key);
void DeleteBst(void);
void DeleteBst(shared_ptr<BstType> Node);
void PrintBst(void);
void PrintBst(const shared_ptr<BstType> Node);
int GetSize(const shared_ptr<BstType> Node);
int GetMin(void);
shared_ptr<BstType> GetMin(const shared_ptr<BstType> Node);
int GetMax(void);
void DeleteMin(void);
shared_ptr<BstType> DeleteMin(shared_ptr<BstType> Node);
int GetRank(const int key);
int GetRank(shared_ptr<BstType> Node, const int key);

//Knights tour (and N-queens problem)
#define CHESS_BOARD_SIZE 8 //chess board size NxN
#define CHESS_RC_TO_BOARD(R,C) ((R*CHESS_BOARD_SIZE)+C)
vector<pair<int, int>> KnightDirections = { { -1,-2 },{ -2,-1 },{ 1,2 },{ 2,1 },
    { -1,2 },{ 2,-1 },{ 1,-2 },{ -2,1 } };
vector<pair<int, int>> PreTraverse, PostTraverse;
bitset<CHESS_BOARD_SIZE*CHESS_BOARD_SIZE> SqVisited;
vector<int> SqParent(CHESS_BOARD_SIZE*CHESS_BOARD_SIZE, -1);
int LastSqVisited = 0;
void KnightsTourAlgm(const int R, const int C);

//N - queens problem
int NqueensLoopCtr = 0;
vector<pair<int, int>> NqueensLocVec;
vector<vector<bool>> NextQueenLocTable (CHESS_BOARD_SIZE, vector<bool>(CHESS_BOARD_SIZE, true));
void NqueensProblem(void);
bool IsValidChessSquare(const int R, const int C);
void UpdateNextQueenLocTable(void);
bool NqueensProblemSol2(const int Col);
bool IsQueenPosValid(const int Row, const int Col);
int Knapsack_Recursive(const int curridx, const int curr_capacity, const vector<int> &knapsack_val, const vector<int> &knapsack_weights);
bool IsUglyNum(const int n);
int StrToNum(const string &In, const int base);
bool IsBoardPntValid(const pair<int,int> NextPnt, const int M, const int N,
                     const vector<vector<bool>> &InBoard,
                     const vector<vector<bool>> &BoardVisited);
void SwapMatElem(int &Val1, int &Val2);
shared_ptr<BstType> SortedArrayToBst(const vector<int> A, const int lo, const int hi);
void ConvNumToBitsCnt(const int num, vector<int> &BitCnt);
int FindArraySum(const vector<int> &A, const int Start, const int End, int &MaxVal);
int FindTotalPaths(const vector<vector<int>> &SumPathsInput, const int M, const int N, const int row, const int col);
shared_ptr<BstType> BuildTreeFromInPreOrder
(const vector<char> &Inorder, const vector<char> &Preorder, int &PreorderIdx, const int lo, const int hi);
int FindInorderIdx(const vector<char> &Inorder, const char val, const int lo, const int hi);
bool IsSplitArraySumPossible(const vector<int> &A, const int idx, const bool IsB1, int B1sum, int B2sum);
long int CustomPower(const int x, const int y);
int MaxProduct(const vector<int> &A, const int idx, const int cnt);
void LogBishopPath(vector<vector<int>> &Temp, const int R, const int C, const int N);
bool PairSortComp(pair<int,int> i, pair<int,int> j);
void IslandsBfs(const int row, const int col, vector<vector<bool>> &Visited, const vector<vector<int>> &IslandsInput);
int LargestNonAdjSum(const vector<int> &InputVec, const int idx, unordered_map<int,int> &Memo);
void AddElemMedian(const int i, priority_queue<int, vector<int>, greater<int>> &MinHeap, priority_queue<int> &MaxHeap);
void BalanceMedian(priority_queue<int, vector<int>, greater<int>> &MinHeap, priority_queue<int> &MaxHeap);
float GetMedian(const priority_queue<int, vector<int>, greater<int>> &MinHeap, const priority_queue<int> &MaxHeap);
void StringPerm(const vector<string> &Input, const string::size_type idx, const string &input_str, unordered_set<string> &StrPermOutputMap);

int main(int argc, const char * argv[])
{
    // insert code here...
    /*
     std::cout << "Hello, World!" << std::endl;
     std::cout << "Enter a number and a string: ";
     long int l_Test = 33;
     char s_Temp[20];
     std::cin >> l_Test >> s_Temp;
     std::cout << "Entered number and string: " << l_Test << " " << s_Temp << std::endl;
     int k = 0;
     for( int i = 0 ; i < 2 ; ++i, k++ )
        std::cout << i << " " << k << std::endl;
     std::cout << k << std::endl;
     */
    /*
     std::cout << "Enter numbers: " << std::endl;
     unsigned long int q_Input = 0, q_Sum = 0;
     if( std::cin >> q_Input )
     {
         while( std::cin >> q_Input )
         {
             std::cout << "Input number: " << q_Input << std::endl;
             q_Sum += q_Input;
         }
     }
     std::cout << "Sum :" << q_Sum << std::endl;
     */
    /*
     Sales_item book1, book2;
     cout << "Enter book1 and book 2: ";
     cin >> book1 >> book2;
     cout << "Book1 details " << book1 << endl;
     cout << "Book2 details " << book2 << endl;
     cout << "Book sum" << book1 + book2 << endl;
     while( cin >> book1 )
        cout << "Book loop " << book1 << endl;
     cerr << "Test error msg" << endl;
     */
    /*unsigned int q_A = 4, q_B = 10;
     int l_C = -40;
     cout << "A-B " << q_A - q_B << " B-A " << (int)q_B + l_C << endl;
     
     unsigned u = 10, u2 = 42;
     std::cout << u2 - u << std::endl;
     std::cout << u - u2 << std::endl;
     int i = 10, i2 = 42;
     std::cout << i2 - i << std::endl;
     std::cout << i - i2 << std::endl;
     std::cout << i - u << std::endl;
     
     char test_char[] = "A";
     cout << test_char << " test msg" " continue " << sizeof(test_char) << endl;
     */
    /*
     long double ld = 3.1415926536;
     int a(ld), b = (ld); // error: narrowing conversion required
     int c(ld), d = ld; // ok: but value will be truncated
     
     int units_sold1 = 0;
     int units_sold2 = {0};
     int units_sold3{0};
     int units_sold4(0);
     string test_string = "test string sdsfsffsdfsfsfsfsdfsdfsd";
     cout << test_string << " test msg" " continue " << sizeof(test_string) << endl;
     
     double x, y = 10.9;
     x = y = 10.9;
     auto xyz = pow(2,32) - 1;
     decltype(xyz) a2 = -4.7;
     cout << sizeof(xyz) << " " << sizeof(a2) << " " << xyz << " " << a2 << endl;
     
     double d1_img;
     */
    /*
     int unique = 0; // unique has block scope
     // output #1: uses global reused; prints 42 0
     std::cout << reused << " " << unique << std::endl;
     int reused = 0; // new, local object named reused hides global reused
     // output #2: uses local reused; prints 0 0
     std::cout << reused << " " << unique << std::endl;
     // output #3: explicitly requests the global reused; prints 42 0
     std::cout << ::reused << " " << unique << std::endl;
     */
    /*
     int i = 100, sum = 0;
     for (int i = 0; i != 10; ++i)
        sum += i;
     std::cout << i << " " << sum << std::endl;
     */
    /*
     int i = 1024, i2 = 2048; // i and i2 are both ints
     int &r = i, r2 = i2; // r is a reference bound to i; r2 is an int
     int i3 = 1024, &ri = i3; // i3 is an int; ri is a reference bound to i3
     int &r3 = r, &r4 = i2; // both r3 and r4 are references
     r3 = 10;
     cout << i << " " << sizeof(i) << " " << sizeof(r3) << endl;
     */
    /*
     int i = 1024, *p = &i, &r = i;
     int *&r1 = p;
     *r1 = 10;
     cout << i << endl;
     
     int ival = 1024;
     int *pi = &ival; // pi points to an int
     int **ppi = &pi; // ppi points to a pointer to an int
     */
    /*
     int i = 42;
     const int &r1 = i; // we can bind a const int& to a plain int object
     const int &r2 = 48; // ok: r1 is a reference to const
     const int &r3 = r1 * 2; // ok: r3 is a reference to const
     i++;
     cout << i << " " << r1 << " " << r2 << " " << r3 << endl;
     
     int a = 10;
     const int *ptr = &a;
     a++;
     cout << a << " " << *ptr << endl;
     */
    /*
     int i = -1, &r = 0, i2 = 2;
     int *const p2 = &i2;
     const int i5 = -1, &r5 = 0;
     const int *const p3 = &i2;
     const int *p1 = &i2;
     const int &const r2;
     const int i3 = i, &r4 = i;
     */
    /*
     int i = 0;
     int *const p1 = &i; // we can't change the value of p1; const is top-level
     const int ci = 42; // we cannot change ci; const is top-level
     const int *p2 = &ci; // we can change p2; const is low-level
     const int *const p3 = p2; // right-most const is top-level, left-most is not
     const int &r = ci; // const in reference types is always low-level
     
     i = ci; // ok: copying the value of ci; top-level const in ci is ignored
     p2 = p3; // ok: pointed-to type matches; top-level const in p3 is ignored
     
     int *p = p3; // error: p3 has a low-level const but p doesn't
     p2 = p3; // ok: p2 has the same low-level const qualification as p3
     p2 = &i; // ok: we can convert int* to const int*
     const int &r3 = ci; // error: can't bind an ordinary int& to a const int object
     const int &r2 = i; // ok: can bind const int& to plain int
     
     typedef double wages; // wages is a synonym for double
     wages base, *p; // base is a synonym for double, p for double*
     */
    /*
     int i = 0, &r = i;
     const int ci = i, &cr = ci;
     auto k = ci, &l = i; // k is int; l is int&
     auto &m = ci, *p = &ci; // m is a const int&;p is a pointer to const int
     // error: type deduced from i is int; type deduced from &ci is const int
     auto &n = i;
     auto *p2 = &ci;
     */
    /*
     // decltype of an expression can be a reference type
     int i = 42, *p = &i, &r = i;
     decltype(r + 0) b; // ok: addition yields an int; b is an (uninitialized) int
     decltype(*p) c; // error: c is int& and must be initialized
     */
    /*int i = 10;
     int &r = i;
     i++;
     r++;
     cout << i << " " << r << endl;
     
     constexpr const int *ptr = nullptr;
     (*ptr)++;
     */
    /*
     string str(10,'q'), str2;
     cout << str << " " << str.size() << " " << sizeof(str) << endl;
     //cin >> str >> str2;
     //cout << str << " " << str2 << " " << str.size() << " " << sizeof(str) << endl;
     
     //while( cin >> str )
     //cout << str << " " << str.size() << " " << sizeof(str) << endl;
     
     //while( getline(cin, str) )
     getline(cin, str);
     if( !str.empty() )
        cout << str << " " << str.size() << " " << sizeof(str) << endl;
     
     for( auto c : str )
        cout << c << " " << " " << sizeof(c) << endl;
     */
    /*
     string str;
     getline(cin, str);
     cout << str << " " << str.size() << " " << endl;
     string temp;
     for( auto c : str )
     {
         if ( !isspace(c) && isalnum(c) )
         {
            temp += c;
         }
         else
         {
             if( temp.size() )
             {
                cout << temp << " " << temp.size() << " " << endl;
             }
             temp = "";
         }
     }
     cout << temp << " " << temp.size() << " " << endl;
     */
    /*
     string s("some string");
     for (decltype(s.size()) index = 0; index != s.size() && !isspace(s[index]); ++index)
        s[index] = toupper(s[index]); // capitalize the current     character
     cout << s << endl;
     */
    /*
     string hex_num = "0123456789ABCDEF", output;
     string::size_type n;
     while( cin >> n )
     {
         if( n < hex_num.size() )
         output += hex_num[n];
     }
     cout << output << " " << output.size() << endl;
     */
    /*
     string c("test"), c2 = "test";
     string c3{"test"}, c4{'c','r', 'x'};
     cout << c << c2 << c3 << c4 << c4.size() << endl;
     */
    /*
     vector<string> v5{"hi"}; // list initialization: v5 has one element
     //vector<string> v6("hi"); // error: can't construct a vector from a string literal
     vector<string> v7{"1"}; // v7 has ten default-initialized elements
     vector<string> v8(10, "hi"); // v8 has ten elements with value "hi"
     cout << v5.size() << " " << v7.size() << " " << v8.size() << endl;
     cout << v5[0] << " " << v7[0] << " " << v8[0] << endl;
     return 0;
     
     int i = 10;
     auto &a = i;
     */
    /*
     string word;
     vector<string> svec;
     while (cin >> word)
     {
        svec.push_back(word);
     }
     //for( decltype(svec.size()) i = 0; i < svec.size(); ++i )
     for( auto i : svec )
        cout << i << " ";
     cout << endl << svec.size() << endl;
     */
    /*
     vector<int> v1 = {1,2,3,4}, v2 = {1,2,4,3};
     cout << (v1 == v2) << endl;
     */
    /*
     vector<int> input = {42, 65, 95, 100, 39, 67, 95, 76, 88, 76, 83, 92, 76, 93};
     vector<int> result(11,0);
     for( const auto i : input )
        if( i <= 100 )
            ++result[i/10];
     for( const auto i : result )
        cout << i << " ";
     cout << endl;
     */
    /*string s("test world");
     for( auto i = s.begin(); i != s.end() && !isspace(*i); ++i )
        *i = toupper(*i);
     cout << s << endl;
     
     int i = 10, *const ptr = &i;
     *ptr = 11;
     (*ptr)++;
     */
    /*
     vector<string> svec;
     string word;
     while( cin >> word )
        svec.push_back(word);
     for( const auto s : svec )
        //if( !s.empty() )
        if( isalnum( *s.cbegin() ) )
        {
            cout << s << " ";
            break;
        }
     cout << endl;
     */
    /*
     int i = 2 + 4 / 2;
     cout << i << endl;
     */
    /*
     unsigned cnt = 42; // not a constant expression
     constexpr unsigned sz = 42; // constant expression
     // constexpr see § 2.4.4 (p. 66)
     int arr[10]; // array of ten ints
     int *parr[sz]; // array of 42 pointers to int
     string bad[cnt]; // error: cnt is not a constant expression
     
     char a1[] = {'C', '+', '+'}; // list initialization, no null
     char a2[] = {'C', '+', '+', '\0'}; // list initialization, explicit null
     char a3[] = "C++"; // null terminator added     automatically
     //const char a4[6] = "Daniel"; // error: no space for the null!
     
     char ca[] = "C++\0";//{'C', '+', '+'}; // not null terminated
     cout << strlen(ca) << endl; // disaster: ca isn't null terminated
     */
    /*
     int i = 0, a[4][5];
     for( auto &row : a )
        for( auto &col : row )
            col = i++;
     
     for( const auto &row : a )
     {
         for( const auto col : row )
            cout << col << " ";
         cout << endl;
     }
     
     for( auto i = begin(a); i != end(a); ++i )
        for( auto j = begin(*i); j != end(*i); ++j )
            cout << *j << " ";
     cout << endl;
     
     for( auto i = begin(a); i != end(a); ++i )
        cout << *((int *)(*i) + 2) << " ";
     cout << endl;
     
     for( auto i = begin(a[3]); i != end(a[3]); ++i )
        cout << *i << " ";
     cout << endl;
     */
    /*
     int i;
     while((cin >> i) && (99 != i))
        cout << i << " ";
     cout << "exit loop" << endl;
     */
    /*
     typedef int int16;
     int16 a[10];
     using b = int[10];
     b c;
     cout << sizeof(b) << endl;
     */
    /*
     Book b1;
     Book b2("Harry Potter", 15, 10);
     b1.Change("Perry Mason", b2.GetPrice());
     Book b3;
     Book b4(b2);
     
     Read(cin, b3);
     
     Print(cout, b1);
     Print(cout, b2);
     cout << b3.GetTitle() << " " << b3.GetPrice() << " " << b3.GetNumSold() << endl;
     Print(cout, b4);
     
     cout << "End Test Class" << endl;
     */
    /*
     ifstream fin;
     string filename = "/Users/hari/C-programs/Xcode/TestCPP/TestCPP/TestCPP/TestInput.txt";
     fin.open(filename);
     vector<int> ivec;
     int temp;
     
     if(!fin)
        cerr << "Unable to open input file" << endl;
     else
     {
         cout << "file open success" << endl;
         while( fin >> temp )
            ivec.push_back(temp);
         fin.close();
         for(auto i : ivec)
            cout << i << endl;
     }
     
     ofstream fout;
     fout.open(filename + ".copy");
     if(!fout)
        cerr << "Unable to open output file" << endl;
     else
     {
        for(auto i = ivec.begin(); i != ivec.end(); ++i)
        fout << *i << endl;
     }
     fout.close();
     */
    /*
     //vector<int> ivec{2,4,5,6,8,10,14,15,19,20,30,40,50,60};
     vector<int> ivec{87,5,690,87,87,74,2,4,5,87,10,18,765,87,16,5,8,5,3,2,11};
     int search = 18;
     decltype(ivec.size()) e, b, m;
     bool elem_found = false;
     int cnt = 0;
     b = 0;
     e = ivec.size();
     m = (b+e)/2;
     while(b<e)
     {
         ++cnt;
         if(ivec[m] == search)
         {
             cout << "Element " << search << " found at index " << m << ". Cnt " << cnt << endl;
             elem_found = true;
             break;
         }
         else if(search < ivec[m])
            e = m;
         else
            b = m+1;
         m = (b+e)/2;
     }
     
     if( false == elem_found )
        cout << "Element " << search << " not found. Cnt " << cnt << endl;
     */
    /*
     vector<string> v = {"quasi", "simba", "frollo", "scar"};
     list<string> slist;
     // insert the last two elements of v at the beginning of slist
     cout << *slist.insert(slist.begin(), v.end() - 2, v.end()) << endl;
     slist.insert(slist.end(), {"these", "words", "will", "go", "at", "the", "end"});
     // run-time error: iterators denoting the range to copy from
     // must not refer to the same container as the one we are changing
     slist.insert(slist.begin(), slist.begin(), slist.end());
     */
    /*
     list<string> lst;
     string word;
     auto iter = lst.begin();
     while (cin >> word)
        iter = lst.insert(iter, word); // same as calling push_front
     */
    /*
     vector<int> ivec {1,2,3,4,5,6,7,8,9,10};
     auto iter = ivec.begin();
     while (iter != ivec.end())
     {
         if((*iter) % 2)
            iter = ivec.erase(iter);
         else
         {
             iter = ivec.insert(iter, *iter);
             iter += 2;
         }
     }
     */
    /*
     const char *cp = "Hello World!!!"; // null-terminated array
     char noNull[] = {'H', 'i'}; // not null terminated
     string s1(cp); // copy up to the null in cp; s1 == "Hello World!!!"
     string s2(noNull,2); // copy two characters from no_null; s2 == "Hi"
     string s3(noNull); // undefined: noNull not null terminated
     string s4(cp + 6, 5);// copy 5 characters starting at cp[6]; s4 == "World"
     string s5(s1, 6, 5); // copy 5 characters starting at s1[6]; s5 == "World"
     string s6(s1, 6); // copy from s1 [6] to end of s1; s6 == "World!!!"
     string s7(s1,6,20); // ok, copies only to end of s1; s7 == "World!!!"
     string s8(s1, 16); // throws an out_of_range exception
     */
    /*
     string s2 = "pi = 3.14";
     double d;
     // convert the first substring in s that starts with a digit, d = 3.14
     d = stod(s2.substr(s2.find_first_of("+-.0123456789")));
     */
    /*
     string line ("the quick red fox jumps over the slow red turtle back");
     vector<string> words;
     string::size_type b = 0;
     while( b < line.size() )
     {
         string::size_type space = line.find(" ",b);
         if(space != string::npos)
            words.push_back(line.substr(b, space-b));
         else
         {
             words.push_back(line.substr(b, line.size()-b));
             break;
         }
         b = space + 1;
     }*/
    
    /*
     vector<string> sort_words(words);
     sort(sort_words.begin(), sort_words.end());
     vector<string> sort_words_copy(sort_words);
     auto end_unique1 = unique(words.begin(), words.end());
     auto end_unique2 = unique(sort_words.begin(), sort_words.end());
     */
    
    /*
     //vector reverse
     vector<string> reverse_words(words);
     auto first = reverse_words.begin(), last = reverse_words.end();
     while((first!=last) && (first!=--last))
     {
         iter_swap(first,last);
         ++first;
     }
     */
    
    /*
     //string reverse
     string reverse_line(line);
     auto first = reverse_line.begin(), last = reverse_line.end();
     while((first!=last) && (first!=--last))
     {
         iter_swap(first,last);
         ++first;
     }
     */
    /*
     string line ("testline the quick/red fox*jumps over*the/slow red turtle back");
     //string line ("the quick red");
     string reverse_line;
     string::size_type end = line.size() - 1;
     while( (end > 0) && (end < line.size()) )
     {
         string::size_type space = line.find_last_of(" /",end);
         if(space != string::npos)
         {
             reverse_line += line.substr(space+1, end-space);
             //reverse_line += string(" ");
             reverse_line += line[space];
         }
         else
         {
             reverse_line += line.substr(0, end+1);
             break;
         }
         end = space - 1;
     }
     */
    /*
     vector<string> svec {"to", "another"};
     sort(svec.begin(), svec.end());
     */
    /*
     istream_iterator<int> in_iter(cin), eof; // read ints from cin
     vector<int> vec(in_iter, eof); // construct vec from an iterator range
     */
    
    /*
     // count the number of times each word occurs in the input
     map<string, int> word_count; // empty map from string to size_t
     string word;
     while (cin >> word)
        ++word_count[word]; // fetch and increment the counter for word
     for (const auto &w : word_count) // for each element in the map
        // print the results
        cout << w.first << " occurs " << w.second << ((w.second > 1) ? " times" : " time") << endl;
     */
    /*
     // count the number of times each word occurs in the input
     map<string, size_t> word_count; // empty map from string to size_t
     set<string> exclude = {"The", "But", "And", "Or", "An", "A",
                             "the", "but", "and", "or", "an",
                             "a"};
     string word;
     while (cin >> word)
     // count only words that are not in exclude
        if (exclude.find(word) == exclude.end())
            ++word_count[word]; // fetch and increment the counter for word
     
     for (const auto &w : word_count) // for each element in the map
        // print the results
        cout << w.first << " occurs " << w.second << ((w.second > 1) ? " times" : " time") << endl;
     */
    /*
     // define a vector with 20 elements, holding two copies of each number from 0 to 9
     vector<int> ivec;
     for (vector<int>::size_type i = 10; i != 0; --i) {
        ivec.push_back(static_cast<int>(i));
        ivec.push_back(static_cast<int>(i)); // duplicate copies of each number
     }
     // iset holds unique elements from ivec; miset holds all 20 elements
     set<int> iset(ivec.cbegin(), ivec.cend());
     multiset<int> miset(ivec.cbegin(), ivec.cend());
     cout << ivec.size() << endl; // prints 20
     cout << iset.size() << endl; // prints 10
     cout << miset.size() << endl; // prints 20
     */
    /*
     // count the number of times each word occurs in the input
     map<string, int> word_count; // empty map from string to size_t
     string word;
     while (cin >> word)
     {
         auto ret = word_count.insert(make_pair(word,1));
         if( false == ret.second )
         ++ret.first->second;
     }
     for (const auto &w : word_count) // for each element in the map
        // print the results
        cout << w.first << " occurs " << w.second << ((w.second > 1) ? " times" : " time") << endl;
     */
    /*
     map<string, int> word_count; // empty map from string to size_t
     cout << word_count["Anna"]; // fetch the element indexed by Anna; prints     1
     ++word_count["Anna"]; // fetch the element and add 1 to it
     cout << word_count["Anna"]; // fetch the element and print it; prints 2
     auto ret = word_count.find("Anna");
     cout << ret->second << (word_count.find("Anna"))->second;
     */
    /*
     ifstream rules("/Users/hari/C-programs/Xcode/TestCPP/TestCPP/TestCPP/rules.txt");
     ifstream input("/Users/hari/C-programs/Xcode/TestCPP/TestCPP/TestCPP/input.txt");
     ofstream output("/Users/hari/C-programs/Xcode/TestCPP/TestCPP/TestCPP/output.txt");
     
     //generate rules map
     map<string, string> rules_map;
     #if 0
     string rules_line;
     while(getline(rules, rules_line))
     {
         istringstream stream_from_line(rules_line);
         string key, value;
         while((stream_from_line >> key) && (stream_from_line >> value))
            rules_map[key] = value;
     }
     #endif
     
     string key, value;
     while((rules >> key) && (rules >> value))
        rules_map[key] = value;
     
     //run through input file
     string input_line;
     while(getline(input, input_line))
     {
         istringstream stream_from_line(input_line);
         string word, new_word;
         while(stream_from_line >> word)
         {
             new_word = word;
             auto ret = rules_map.find(word);
             if(ret != rules_map.end())
                new_word = ret->second;
             output << new_word << " ";
         }
         output << endl;
     }
     
     rules.close();
     input.close();
     output.close();
     */
    /*
     // count the number of times each word occurs in the input
     unordered_map<string, int> word_count; // empty map from string to size_t
     string word;
     while (cin >> word)
        ++word_count[word]; // fetch and increment the counter for word
     for (const auto &w : word_count) // for each element in the map
        // print the results
        cout << w.first << " occurs " << w.second << ((w.second > 1) ? " times" : " time") << endl;
     */
    /*
     auto p = make_shared<ifstream>("/Users/hari/C-programs/Xcode/TestCPP/TestCPP/TestCPP/rules.txt");
     p->close();
     */
    /*
     int *pi = new int(1024); // object to which pi points has value 1024
     string *ps = new string(10, '9'); // *ps is "9999999999"
     // vector with ten elements with values from 0 to 9
     vector<int> *pv = new vector<int>{0,1,2,3,4,5,6,7,8,9};
     */
    /*
     string word;
     vector<string> svec;
     
     cout << "enter words. type 'end' to terminate" << endl;
     while (cin>>word)
     {
         if (word == "end")
            break;
         svec.push_back(word);
     }
     
     if( !svec.empty() )
        for( auto &c : svec )
            c[0] = toupper(c[0]);
     
     for( const auto &s : svec )
        cout << s << " ";
     cout << endl;
     */
    /*
     const int sum = 8;
     vector<int> input_vec {1,2,3,5,8,9,10};
     cout << "enter num. type '-999' to terminate" << endl;
     
     int num;
     while(cin >> num)
     {
         if(-999 == num)
         break;
         input_vec.push_back(num);
     }
     */
    /*
     for (const auto i : input_vec)
        if( true == BinSearch(sum-i, &input_vec))
            cout << "element" << i << endl;
     */
    /*
     int ia[4][4] = { {0,1,2,3}, {4,5,6,7}, {8,9,10,11}, {12,13,14,15} };
     for( auto i:ia[2])
        cout << i << " ";
     cout << endl;
     for(auto i:ia)
        cout << *i << " ";
     cout << endl;
     for(auto i:ia)
        cout << *((int *)i + 2) << " ";
     cout << endl;
     for(int i = 0; i < 4; ++i)
        cout << **(ia + i) << " ";
     cout << endl;
     
     vector<int> ivec(begin(ia[0]), end(ia[0]));
     for(auto i:ivec)
        cout << i << " ";
     cout << endl;
     */
    /*
     cout << g() << f()+ (g()*h()) + j() << j() << endl;
     int i = 0;
     cout << i << " " << i++ << endl;
     
     short short_value = 32767;
     short_value += 1;
     cout << "short_value: " << short_value << endl;
     */
    /*
     ifstream inputf("/Users/hari/C-programs/Xcode/TestCPP/TestCPP/TestCPP/input_para.txt");
     ofstream outputf("/Users/hari/C-programs/Xcode/TestCPP/TestCPP/TestCPP/output_para_srch_results.txt");
     
     string input_line, input_word;
     int linenum = 0, searchcnt = 0, wordnum = 0;
     //const string search_word = "memory";
     set<string> search_word = {"memory", "objects"};
     queue<string> WordQ;
     
     while (getline(inputf,input_line))
     {
     ++linenum;
     wordnum = 0;
     istringstream inputs_stream(input_line);
     while(inputs_stream >> input_word)
     {
         ++wordnum;
         WriteToQueue(WordQ, input_word);
         //trim input word to remove non-alphanumberic characters
         auto non_alnum = input_word.find_first_of(",.:;");
         if( string::npos != non_alnum )
            input_word.resize(non_alnum);
         if( search_word.end() != search_word.find(input_word) )
         {
            ++searchcnt;
            outputf << "  line(" << linenum << ") " << "word(" << wordnum << ") ";
            while(!WordQ.empty())
            {
                outputf << WordQ.front() << " ";
            WordQ.pop();
            }
            outputf << endl;
            }
         }
     }
     outputf << *search_word.cbegin() << "/" << *search_word.crbegin() << " occurs " << searchcnt << " times." << endl;
     
     inputf.close();
     outputf.close();
     */
    
    /*cout << RC_TO_INT(2,4) << endl;
     string sint("a1s2d4");
     cout << stoi(sint)+1 << endl;
     */
    /*
     classA c1(4), c2(c1), c4;
     c1.setTemp(9);
     classA c3 = c1 + c2;
     c3 += c1;
     c4 = c3;
     bool a = c1==c2;
     bool b = c1!=c2;
     bool c = c1<c2;
     cout << a << b << c << endl;
     */
    /*
     unordered_map<int,shared_ptr<classA>> class_map;
     int temp;
     while( cin >> temp )
     {
         auto temp_sptr = make_shared<classA>(temp);
         class_map[temp] = temp_sptr;
     }
     
     //for(const auto &i:class_map)
     //cout << i.first << " " << i.second << endl;
     for(auto i = class_map.cbegin(); i != class_map.cend(); ++i)
        cout << (*i).first << " " << (*i).second << endl;
     */
    /*
     //GCD
     int p,q,r=0,loopcnt=0;
     cin >> p >> q;
     while(0 != q)
     {
         ++loopcnt;
         r = p % q;
         p = q;
         q = r;
     }
     cout << "gcd: " << p << " loopcnt: " << loopcnt << endl;
     */
    /*
     vector<int> A;
     const int S = 88;
     const int N = 100;
     for(auto i = 0; i < N; ++i)
        A.push_back(2*i);
     auto hi = A.size();
     decltype(hi) lo = 0;
     cout << BinSearchRecursive(S, A, lo, hi) << endl;
     */
    /*
     int N;
     cin >> N;
     cout << Factorial(N) << endl;
     */
    
    /*
     //Djikstra's algm to solve expression with braces
     string expr = "(((2+(3*4))/2)";
     stack<int> oprand;
     stack<char> oprtr;
     
     for(const auto &i : expr)
     {
         if(isdigit(i))
            oprand.push(stoi(&i));
         else if(('*'==i) || ('/'==i) || ('+'==i) || ('-'==i))
            oprtr.push(i);
         else if(')'==i)
         {
             int a, b;
             char c;
             if(!Popstack(oprand,a))
             {
                 cout << "Operand stack empty" << endl;
                 break;
             }
             if(!Popstack(oprand,b))
             {
                 cout << "Operand stack empty" << endl;
                 break;
             }
             if(!Popstack(oprtr,c))
             {
                 cout << "Operator stack empty" << endl;
                 break;
             }
             switch(c)
             {
                 case '+':
                 {
                     oprand.push(a+b);
                     break;
                 }
                 case '-':
                 {
                     oprand.push(b-a);
                     break;
                 }
                 case '*':
                 {
                     oprand.push(a*b);
                     break;
                 }
                 case '/':
                 {
                     oprand.push(b/a);
                     break;
                 }
             }
         }
         else if('('==i)
         {
            //ignore
         }
         else
         {
            cout << "Invalid character: " << i << endl;
            break;
         }
     }
     
     int sol;
     bool ret = Popstack(oprand,sol);
     cout << "Solution: " << sol << " " << ret << endl;
     */
    
    /*
     //2SUM
     vector<int> ivec = {2,6,4,9,6,3,5,2,2,2,14,-3};
     int sum = 11;
     unordered_multiset<int> iset;
     
     for(const auto &i : ivec)
     {
         iset.insert(i);
         if(iset.end()!=iset.find(sum-i))
         {
             auto myrange = iset.equal_range(sum-i);
             while(myrange.first != myrange.second)
             {
                 cout << "Pair found: " << i << " " << sum-i << endl;
                 ++myrange.first;
             }
         }
     }
     */
    
    /*
     //3SUM
     vector<int> ivec = {2,6,4,9,6,3,5,2,2,2,14,-3};
     int sum = 11;
     unordered_multiset<int> iset;
     
     for(auto i = ivec.cbegin(); i != ivec.cend(); ++i)
     {
         iset.insert(*i);
         for(auto j = i+1; j < ivec.cend(); ++j)
         {
             int find_elem = sum-*i-*j;
             if(iset.end()!=iset.find(find_elem))
             {
                 auto myrange = iset.equal_range(find_elem);
                 while(myrange.first != myrange.second)
                 {
                     cout << "3SUM Pair found: " << *i << " " << *j << " " << find_elem << endl;
                     ++myrange.first;
                 }
             }
         }
     }
     */
    
    /*
     //Union find
     list<unordered_set<int> > setlist;
     const int N = 20;
     cout << "Max nodes: 1 to " << N << endl;
     cout << "1 Node#1 Node#2 for connecting" << endl;
     cout << "2 Node#1 Node#2 for checking" << endl;
     cout << "3 Print" << endl;
     int op = 0, node1 = 0, node2 = 0;
     while((cin >> op) && ((3 == op) || ((cin >> node1) && (cin >> node2))))
     {
         if((node2>N)||(node2>N)||(op<1)||(op>3))
         {
             cerr << "Invalid op: " << op << " Nodes: " << node1 << ", " << node2 << endl;
             break;
         }
     
         if(1==op)
         {
             auto node1itr = setlist.end(), node2itr = setlist.end();
             for(auto i = setlist.begin(); i != setlist.end(); ++i)
             {
                 const unordered_set<int> &tempset = *i;
     
                 if(tempset.end() != tempset.find(node1))
                    node1itr = i;
                 if(tempset.end() != tempset.find(node2))
                    node2itr = i;
                 if((node1itr != setlist.end()) && (node2itr != setlist.end()))
                    break;
             }
             if((node1itr != setlist.end()) && (node2itr != setlist.end()))
             {
                 unordered_set<int> &tempset1 = *node1itr, &tempset2 = *node2itr;
     
                 // Merge node2 set onto node1 set. Delete node2 set
                 if(node1itr != node2itr)
                 {
                     tempset1.insert(tempset2.cbegin(), tempset2.cend());
                     setlist.erase(node2itr);
                 }
                 else
                    cout << node1 << " and " << node2 << " are already connected" << endl;
             }
             else if((node1itr == setlist.cend()) && (node2itr != setlist.cend()))
             {
                 unordered_set<int> &tempset2 = *node2itr;
                 tempset2.insert(node1); //Add node1 to node2 set
             }
             else if((node1itr != setlist.cend()) && (node2itr == setlist.cend()))
             {
                 unordered_set<int> &tempset1 = *node1itr;
                 tempset1.insert(node2); //Add node2 to node1 set
             }
             else
             {
                 //Node1 and node2 not present in setlist. Create new set and add to list
                 unordered_set<int> tempset = {node1,node2};
                 setlist.push_back(tempset);
             }
         }
         else if(2==op)
         {
             auto i = setlist.begin();
             for(; i != setlist.end(); ++i)
             {
                 const unordered_set<int> &tempset = *i;
                 if(tempset.end() != tempset.find(node1))
                 {
                     //Found set with node1. Check if node2 is present in this set, and return
                     if(tempset.end() != tempset.find(node2))
                        cout << node1 << " and " << node2 << " are connected" << endl;
                     else
                        cout << node1 << " and " << node2 << " are NOT connected" << endl;
                     break;
                 }
                 if(tempset.end() != tempset.find(node2))
                 {
                     //Found set with node2. Check if node1 is present in this set, and return
                     if(tempset.end() != tempset.find(node1))
                        cout << node1 << " and " << node2 << " are connected" << endl;
                     else
                        cout << node1 << " and " << node2 << " are NOT connected" << endl;
                     break;
                 }
             }
             if(setlist.end() == i)
             {
                //Reached end of 'i' loop. This means Node1 and node2 were not found in the setlist
                cout << node1 << " and " << node2 << " are NOT connected" << endl;
             }
         }
         else if(3==op)
         {
            //Print setlist
            PrintSetlist(setlist);
         }
     
         //Clear op for while loop
         op = 0;
     }
     */
    
    /*
     //Selection sort
     //Find the min. element in the array and place it at the beginning
     //Const complexity of O(N) even if input is already sorted
     //typedef int SortType;
     //vector<SortType> ivec = {3,16,9,7,1,7,7,87,45,76};
     //vector<SortType> ivec = {"cat", "and", "zoo", "yes", "comic", "yanni", "able", "axle", "history"};
     //SortType temp;
     vector<int> ivec = {3,16,9,7,1,7,7,87,45,76};
     int temp;
     for(auto i = ivec.begin(); i != ivec.end(); ++i)
     {
         //Find min. in leftover array
         auto min = i;
         for(auto j = i; j != ivec.end(); ++j)
            if(*j < *min)
                min = j;
         //Swap minimum with LHS of array
         temp = *min;
         *min = *i;
         *i = temp;
     }
     
     for(const auto &i : ivec)
        cout << i << " ";
     cout << endl;
     */
    
    /*
     //Insertion sort
     //Traverse through each element and place it in the correct postion in the array
     //before it. If vector is already sorted, complexity is linear. No need to run 2nd loop
     //vector<int> ivec = {3,16,9,7,1,7,7,87,45,76};
     vector<int> ivec = {1,3,5,6,7,10,16,18,20,27,29};
     int temp;
     for(auto i = ivec.begin()+1; i != ivec.end(); ++i)
     {
         for(auto j = i; j != ivec.begin() && (*j < *(j-1)); --j)
         {
             temp = *j;
             *j = *(j-1);
             *(j-1) = temp;
             PrintSortVec(ivec);
         }
     }
     */
    
    /*
     //Merge sort
     //Sort 2 halves separately and then combine them. NlogN complexity
     //vector<int> ivec = {1,3,5,6,7,10,87,2,4,5,7,9,18,20,34};
     vector<int> ivec = {3,16,9,7,1,7,7,87,45,76};
     auto hi = ivec.size()-1;
     decltype(hi) lo = 0;
     MyMergeSort(ivec, lo, hi);
     */
    
    /*
     //Quick sort
     //Split/partition the array based on a partitioning element, and then move all items
     //less than the element to the left, and greater than the element to its right.
     //Sort the 2 partitions individually
     vector<int> ivec = {8,16,9,7,1,7,7,87,45,76};
     auto hi = ivec.size()-1;
     decltype(hi) lo = 0;
     MyQuickSort(ivec, lo, hi);
     */
    
    /*
     //Shuffle vector
     vector<int> ivec;
     for(int i = 0; i < 20; ++i)
     {
        ivec.push_back(i);
     }
     for(int i = 0; i < 10; ++i)
     {
         ShuffleVec(ivec);
         PrintSortVec(ivec);
     }
     */
    /*
     vector<int> ivec = {4,2,7,9,10,14,6,66,45,87,23,12,11,1,11,17,3};
     priority_queue<int> pq (ivec.cbegin(), ivec.cend());
     while(!pq.empty())
     {
         cout << pq.top() << " ";
         pq.pop();
     }
     cout << endl;
     */
    
    /*
     //Binary heap tree (max priority tree)
     vector<int> ivec = {4,2,7,9,10,14,6,66,45,87,23,12,11,1,11,17,3};
     ShuffleVec(ivec);
     PrintSortVec(ivec);
     vector<int> heaptree;
     for(const auto &i : ivec)
     {
         InsertHeapTree(heaptree, i);
         PrintSortVec(heaptree);
     }
     while(!heaptree.empty())
     {
         cout << "GetMax: " << GetHeapTree(heaptree) << " " << endl;
         PrintSortVec(heaptree);
     }
     */
    
    /*
     //Heapsort
     vector<int> ivec = {4,2,7,9,10,14,6,66,45,87,23,12,11,1,11,17,3};
     ShuffleVec(ivec);
     PrintSortVec(ivec);
     //construct heap upto half of the size
     for(auto i = ivec.size()/2; i < ivec.size(); --i)
        HeapSwimDown(ivec,i,ivec.size());
     PrintSortVec(ivec);
     
     //Swap Max element with end of array. Reduce array size and construct heap again
     auto vecsize = ivec.size();
     while(vecsize)
     {
         SwapIvec(ivec,0,vecsize-1);
         --vecsize;
         HeapSwimDown(ivec,0,vecsize);
     }
     PrintSortVec(ivec);
     */
    
    /*
     //Search for permutations of a string within another string
     const string s = "";
     const string y = "abcabckdakckkjbkbcaaabcbbc";
     unordered_multiset<char> srchset, srchset_copy;
     for(const auto &i : s)
        srchset.insert(i);
     srchset_copy = srchset;
     int srchcnt = 0;
     for(string::size_type i = 0; i < y.size(); ++i)
     {
         if(srchset_copy.find(y[i]) != srchset_copy.end())
         {
             ++srchcnt;
             auto range = srchset_copy.equal_range(y[i]);
             srchset_copy.erase(range.first); //Erase only 1 instance of character
         }
         else
         {
             srchcnt = 0;
             if(srchset_copy.size() != srchset.size())
             //Make a copy only if its not already the same
                srchset_copy = srchset;
         }
         if(srchcnt == s.size() && !s.empty())
         {
             cout << "Found: " << i << " " << y.substr(i+1-s.size(), s.size()) << endl;
             srchcnt = 0;
             srchset_copy = srchset;
             i += 1;
             i -= s.size(); //Move the loop counter backwards by s.size()-1 to find nested strings
         }
     }
     */
    
    /*
     //MxN connected component matrix - solution #1
     #define GET_NUM(R,C)    (((R)*N)+(C))
     const int M = 7, N = 6;
     const int In[M][N] =   {{0,1,1,0,0,0},
                             {0,0,0,1,1,1},
                             {1,1,1,0,0,0},
                             {0,0,0,1,0,0},
                             {0,1,0,0,0,1},
                             {0,1,0,1,0,0},
                             {0,1,1,1,0,0}};
     list<unordered_set<int>> setlist;
     vector<pair<int,int>> Ngbrs {{-1,0},{0,-1},{1,0},{0,1}};
     
     //Add to setlist
     for(int i = 0; i < M; ++i)
        for(int j = 0; j < N; ++j)
             if(In[i][j])
             {
                 const int CurrNode = GET_NUM(i,j);
                 vector<decltype(setlist.begin())> NgbrItr;
                 for(auto &k : Ngbrs)
                 {
                    int NgbrRow = i+k.first, NgbrCol = j+k.second;
                    if((NgbrRow<0) || (NgbrRow>=M)||
                        (NgbrCol<0) || (NgbrCol>=N))
                        continue;
                    int Ngbr = GET_NUM(NgbrRow, NgbrCol);
                    for(auto l = setlist.begin(); l != setlist.end(); ++l)
                    {
                        unordered_set<int> &TempSet = *l;
                        if(TempSet.find(Ngbr) != TempSet.end())
                            NgbrItr.push_back(l);
                    }
                 }
     
                 if(NgbrItr.empty())
                 {
                    //Input not found. Create a new set (component) and add to the list
                    unordered_set<int> TempSet = {CurrNode};
                    setlist.push_back(TempSet);
                    continue;
                 }
                 if(NgbrItr.size() > 1)
                 {
                     //More than 1 ngbr is already present in the setlist. Merge all of them into 1 set and
                     //delete the set that was merged
                     for(auto l = 1; l < NgbrItr.size(); ++l)
                     {
                         NgbrItr[0]->insert(NgbrItr[l]->begin(), NgbrItr[l]->end());
                         setlist.erase(NgbrItr[l]);
                     }
                 }
     
                 //Add the new element to the setlist
                 NgbrItr[0]->insert(CurrNode);
            }
     */
    
    /*
     //MxN connected component matrix - solution #2
     int ComponentCnt = 0;
     for(int i = 0; i < M; ++i)
        for(int j = 0; j < N; ++j)
             if((ComponentInputs[i][j]) && (-1 == ComponentRoot[i][j]))
             {
                 ++ComponentCnt;
                 SetRoot(i,j, ComponentCnt);
             }
     
     for(int i = 0; i < M; ++i)
     {
         for(int j = 0; j < N; ++j)
            cout << ComponentRoot[i][j] << " ";
         cout << endl;
     }
     */
    
    /*
     //Reverse vowels of a string
     string s = "xcvtr";
     vector<pair<char,int>> VowelPos;
     const unordered_set<char> Vowels = {'a','e','i','o','u'};
     for(string::size_type i = 0; i < s.size(); ++i)
        if(Vowels.find(s[i]) != Vowels.end())
            VowelPos.push_back(make_pair(s[i], i));
     for(decltype(VowelPos.size()) i = 0; i < VowelPos.size(); ++i)
        s[VowelPos[i].second] = VowelPos[VowelPos.size()-1-i].first;
     */
    
    /*
     //Atom counter
     const string formula = "(H2O2He3Mg4O4)2";
     map<string,int> AtomMap;
     string AtomName;
     int AtomCnt = 0;
     bool IsParanthesis = false;
     for(const auto &i : formula)
     {
         if(isalpha(i))
         {
             if(isupper(i))
             {
                 //Start New atom. Before that, store previous atom in map
                 if(!AtomName.empty())
                 {
                     auto MapItr = AtomMap.find(AtomName);
                     if(AtomMap.find(AtomName) == AtomMap.end())
                         //Atom not present. Create new element
                         AtomMap[AtomName] = AtomCnt;
                     else
                         //Atom already exists. Add to previous cnt
                         MapItr->second += AtomCnt;
                 }
                 AtomName = i;
                 AtomCnt = 1;
             }
             else
                 //Atom name has multiple characters
                 AtomName += i;
         }
         if(isdigit(i))
         {
            int Num = static_cast<int>(i) - static_cast<int>('0');
            if(IsParanthesis)
            {
                //Multiply all atom counts in map
                for(auto &j : AtomMap)
                    j.second *= Num;
                IsParanthesis = false;
            }
            else
                AtomCnt = Num;
         }
         if(i == ')')
         {
             IsParanthesis = true;
             //Store previous atom in map
             if(!AtomName.empty())
             {
                 auto MapItr = AtomMap.find(AtomName);
                 if(AtomMap.find(AtomName) == AtomMap.end())
                     //Atom not present. Create new element
                     AtomMap[AtomName] = AtomCnt;
                 else
                     //Atom already exists. Add to previous cnt
                     MapItr->second += AtomCnt;
             }
         }
     }
     */
    
    /*
     //Staircase problem (coin problem)
     const int N = 7; //need to reach N steps using any combination of below steps
     vector<int> Steps = {5,3,1};
     unordered_map<int,int> ResultMap; //Temporary space to store results to avoid multiple recursion
     cout << NumWaysStair(N,ResultMap,Steps) << endl;
     for(const auto &i : ResultMap)
        cout << i.first << " " << i.second << endl;
     */
    
    /*
     //Find the Longest Increasing Subsequence
     //const vector<int> InputVec = {10,8,3,9,2,11,4,7,12,13,14};
     const vector<int> InputVec = {6,4,5,3,2,1,7,8,9,10,11};
     list<set<int>> SetList;
     for(const auto &i : InputVec)
     {
         bool InputAdded = false;
         set<int>::size_type PrevSetSize = 0;
         for(auto j = SetList.begin(); j != SetList.end(); ++j)
         {
             set<int> &TempSet = *j;
             if(i >= *TempSet.rbegin())
             {
                 if(TempSet.size() == PrevSetSize)
                     //Delete this set as the same element is added to 2 sets of same size
                     SetList.erase(j);
                 else
                 {
                     TempSet.insert(i);
                     InputAdded = true;
                     PrevSetSize = TempSet.size() - 1;
                 }
             }
         }
         if(!InputAdded)
            //Create a new set
            SetList.emplace_back(set<int>{i});
     }
     set<int>::size_type MaxSize = 0;
     for(const auto &i : SetList)
     {
         const set<int> &TempSet = i;
         if(MaxSize < TempSet.size())
            MaxSize = TempSet.size();
     }
     for(const auto &i : SetList)
     {
         const set<int> &TempSet = i;
         if(MaxSize == TempSet.size())
         {
             for(const auto &j : TempSet)
                cout << j << " ";
             cout << endl;
         }
     }
     */
    
    /*
     //Num ways to get sum (N-sum)
     const int Sum = 16; //need to reach N steps using any combination of below steps
     vector<int> InputVec = {6,4,10,2};
     unordered_map<string,int> ResultMap; //Temporary space to store results to avoid multiple recursion
     cout << NumWaysNsum(Sum,InputVec.size()-1,ResultMap,InputVec) << endl;
     for(const auto &i : ResultMap)
        cout << i.first << " " << i.second << endl;
     */
    
    /*
     //Binary search tree BST
     cout << "1: Insert key and value" << endl;
     cout << "2: Delete key" << endl;
     cout << "3: Get value[key]" << endl;
     cout << "4: Print" << endl;
     cout << "5: Delete BST" << endl;
     cout << "6: Get min. key" << endl;
     cout << "7: Get rank of a key" << endl;
     cout << "0: Exit" << endl;
     int option = 99;
     while (cin >> option)
     {
         if (0 == option)
            break;
         else if (1 == option)
         {
             int key;
             string value;
             cin >> key >> value;
             InsertKey(key, value);
         }
         else if (2 == option)
         {
             int key;
             cin >> key;
             DeleteKey(key);
         }
         else if (3 == option)
         {
             int key;
             cin >> key;
             cout << "Get: Key: " << key << " Value: " << GetValue(key) << endl;
         }
         else if (4 == option)
            PrintBst();
         else if (5 == option)
            DeleteBst();
         else if (6 == option)
            cout << "Get Min: Key: " << GetMin() << endl;
         else if (7 == option)
         {
             int key;
             cin >> key;
             cout << "Get: Rank: " << GetRank(key) << " Key: " << key << endl;
         }
         else
         {
             cout << "Invalid option: " << option << endl;
             break;
         }
     }
     DeleteBst();
     */
    
    /*
     //Graph
     const int N = 15;
     int TestConnect[N][2] = { { 0,1 }, { 0,2 }, { 3,4 }, { 5,8 }, { 4,0 },
                                 { 8,8 }, { 9,2 }, { 6,4 }, { 2,5 }, { 0,5 },
                                 { 3,0 }, { 7,8 }, { 2,7 }, { 9,4 }, { 6,8 } };
     
     for (int i = 0; i < N; ++i)
     {
         MyDfs.AddConnection(TestConnect[i][0], TestConnect[i][1]);
         MyBfs.AddConnection(TestConnect[i][0], TestConnect[i][1]);
     }
     MyDfs.PrintConnections();
     MyDfs.PerformDfs();
     
     MyBfs.PrintConnections();
     MyBfs.PerformBfs();
     */
    
    /*
     //Knight's tour
     pair<int, int> StartingPnt = make_pair(7, 7);
     int StartingLoc = CHESS_RC_TO_BOARD(StartingPnt.first, StartingPnt.second);
     SqParent[StartingLoc] = StartingLoc;
     KnightsTourAlgm(StartingPnt.first, StartingPnt.second);
     for (int i = 0; i < CHESS_BOARD_SIZE*CHESS_BOARD_SIZE; ++i)
     {
         if (0 == i%CHESS_BOARD_SIZE && i > 0)
         cout << endl;
         cout << setfill('0') << setw(2) << i << " ";
     }
     cout << endl << endl;
     
     cout << "AllSqVisited: " << SqVisited.all() << endl;
     cout << "PreTraverse: ";
     for (const auto &i : PreTraverse)
        cout << CHESS_RC_TO_BOARD(i.first, i.second) << " ";
     cout << endl << endl;
     
     cout << "PostTraverse: ";
     for (const auto &i : PostTraverse)
        cout << CHESS_RC_TO_BOARD(i.first, i.second) << " ";
     cout << endl << endl;
     
     cout << "Path: ";
     int KnightPathSq = LastSqVisited;
     while (KnightPathSq != StartingLoc)
     {
         cout << KnightPathSq << " ";
         KnightPathSq = SqParent[KnightPathSq];
     }
     cout << KnightPathSq << endl << endl;
     */
    
    /*
     //N-queens problem
     NqueensProblemSol2(0);
     
     for (int i = 0; i < CHESS_BOARD_SIZE*CHESS_BOARD_SIZE; ++i)
     {
         if (0 == i%CHESS_BOARD_SIZE && i > 0)
            cout << endl;
         cout << setfill('0') << setw(2) << i << " ";
     }
     cout << endl << endl;
     
     cout << "LoopCtr: " << NqueensLoopCtr << " QueenPos: " << endl;
     for (int i = 0; i < CHESS_BOARD_SIZE; ++i)
     {
         for (int j = 0; j < CHESS_BOARD_SIZE; ++j)
         {
             bool KeyFound = false;
             for (const auto &k : NqueensLocVec)
                 if (k.first == i && k.second == j)
                 {
                     cout << "QQ ";
                     KeyFound = true;
                     break;
                 }
             if(!KeyFound)
                cout << "-- ";
         }
         cout << endl;
     }
     cout << endl;
     */
    
    /*
     //Find largest set of ngbr elements in MxN grid with same color (similar to flood fill)
     #define COLOR_ARRAY_M    4
     #define COLOR_ARRAY_N    4
     typedef pair<int, int> xycoord;
     class ColorArrayType
     {
         public:
             int color;
             set<xycoord> coord;
     };
     vector<ColorArrayType> ColorArrayOutput;
     const vector<vector<int>> ColorInput = {{1,1,2,3},
                                             {1,2,2,2},
                                             {1,1,2,3},
                                             {2,2,1,3}}; //each color represented by a number
     vector<vector<bool>> ColorVisited(COLOR_ARRAY_M, vector<bool>(COLOR_ARRAY_N, false));
     vector<xycoord> NextCoord;
     const vector<xycoord> NgbrCoord = { {1,0},{-1,0},{0,1},{0,-1} };
     
     //Perform DFS for each node, for that particular color
     for(auto i = 0; i != COLOR_ARRAY_M; ++i)
        for (auto j = 0; j != COLOR_ARRAY_N; ++j)
        {
             if (!ColorVisited[i][j])
             {
                 NextCoord.clear(); //vector (Stack) should be empty anyway
                 NextCoord.push_back(make_pair(i, j)); //add first element to stack
     
                 ColorArrayType Temp; //create a new entry in output structure
                 Temp.color = ColorInput[i][j];
     
                 while (!NextCoord.empty())
                 {
                     xycoord CurrCoord = NextCoord.back();
                     NextCoord.pop_back();
                     ColorVisited[CurrCoord.first][CurrCoord.second] = true;
                     Temp.coord.insert(CurrCoord);
     
                     //Add ngbrs that have the same color as the current cell's color. Do so only if its not already
                     //visited yet
                     for (const auto &k : NgbrCoord)
                     {
                         xycoord Ngbr = make_pair(CurrCoord.first + k.first, CurrCoord.second + k.second);
                         if((Ngbr.first >= 0 ) && (Ngbr.first < COLOR_ARRAY_M ) &&
                            (Ngbr.second >= 0) && (Ngbr.second < COLOR_ARRAY_N))
                             if(ColorInput[i][j] == ColorInput[Ngbr.first][Ngbr.second])
                                if(!ColorVisited[Ngbr.first][Ngbr.second])
                                    NextCoord.push_back(Ngbr);
                     }
                 }
     
                ColorArrayOutput.push_back(Temp); //Push the new vector element with connected coordinates to output vector
            }
        }
     */
    
    /*
    //Find max square
    #define MAX_SQUARE_M    5
    #define MAX_SQUARE_N    6
    const vector<vector<int>> MaxSquareInput = { { 1,1,0,0,1,1 },
                                                 { 1,0,1,1,1,1 },
                                                 { 1,1,1,1,1,1 },
                                                 { 1,1,1,1,1,1 },
                                                 { 1,1,1,1,1,1 } };
    vector<vector<int>> MaxSquareWorking (MAX_SQUARE_M, vector<int>(MAX_SQUARE_N,0));
    //Identify end nodes of possible squares
    //Copy 1st row and 1st column to working matrix
    int MaxSquareSize = 0, MaxSquarei = 0, MaxSquarej = 0;
    for (auto i = 0; i != MAX_SQUARE_M; ++i)
        MaxSquareWorking[i][0] = MaxSquareInput[i][0];
    for (auto j = 0; j != MAX_SQUARE_N; ++j)
        MaxSquareWorking[0][j] = MaxSquareInput[0][j];
    for (auto i = 1; i != MAX_SQUARE_M; ++i)
    {
        for (auto j = 1; j != MAX_SQUARE_N; ++j)
        {
            if(!MaxSquareInput[i][j])
                MaxSquareWorking[i][j] = 0;
            else
            {
                MaxSquareWorking[i][j] = min(min(MaxSquareWorking[i][j-1], MaxSquareWorking[i-1][j]),MaxSquareWorking[i-1][j-1]) + 1;
                if( MaxSquareSize <= MaxSquareWorking[i][j] )
                {
                    MaxSquareSize = MaxSquareWorking[i][j];
                    MaxSquarei = i;
                    MaxSquarej = j;
                }
            }
            cout << MaxSquareWorking[i][j] << " ";
        }
        cout << endl;
    }
    */
     
    /*
    //XOR doubly linked list (DLL) prev and next ptrs
    class DLLnodeType
    {
    public:
        int NodeAddr;
        int BothPtr;
    };
    vector<DLLnodeType> DLLnodes;
    const int DLLsize = 8;
    for(int i = 0; i < DLLsize; ++i)
    {
        DLLnodeType Temp;
        Temp.NodeAddr = i+2;
        DLLnodes.push_back(Temp);
    }
    for(auto i = 0; i < DLLnodes.size(); ++i)
    {
        if(0==i)
            DLLnodes[i].BothPtr = 0 ^ DLLnodes[i+1].NodeAddr;
        else if(DLLnodes.size()-1==i)
            DLLnodes[i].BothPtr = 0 ^ DLLnodes[i-1].NodeAddr;
        else
            DLLnodes[i].BothPtr = DLLnodes[i+1].NodeAddr ^ DLLnodes[i-1].NodeAddr;
    }
    int HeadPtr = DLLnodes[0].NodeAddr;
    //Traverse
    int Prev = 0, Curr = HeadPtr, Next = 0;
    for(auto i = 0; i < DLLnodes.size(); ++i)
    {
        Curr = DLLnodes[i].NodeAddr;
        Next = DLLnodes[i].BothPtr ^ Prev;
        Prev = DLLnodes[i].NodeAddr;
        cout << "Curr: " << Curr << " Next: " << Next << " Prev: " << Prev << endl;
    }
    */
    
    /*
    //Autocomplete(AC) from dictionary
    const string AC_input = "d";
    const vector<string> AC_dictionary = {"dust", "elephant", "action", "define", "dear", "deal"};
    //build map (hash table) with all combinations of dictionary
    unordered_multimap<string, string> AC_dictionary_map;
    for (const auto &i : AC_dictionary)
        for (string::size_type j = 0; j < i.size(); ++j)
            AC_dictionary_map.insert(make_pair(i.substr(0, j + 1), i));
    
    auto AC_output = AC_dictionary_map.equal_range(AC_input);
    if(AC_output.first != AC_dictionary_map.end())
        for (auto i = AC_output.first; i != AC_output.second; ++i)
            cout << i->second << endl;
    */
    
    /*
    //Knapsack problem
    const int KNAPSACK_SIZE = 50;
    const int KNAPSACK_CAPACITY = 850;
    //const vector<int> knapsack_val = {5,3,1,3,2};
    //const vector<int> knapsack_weights = {1,2,4,2,5};
    //const vector<int> knapsack_val = {2,3,4,6};
    //const vector<int> knapsack_weights = {2,3,4,5};
    //const vector<int> knapsack_val = {565, 406, 194, 130, 435, 367, 230, 315, 393,
    //125, 670, 892, 600, 293, 712, 147, 421, 255};
    //const vector<int> knapsack_weights = knapsack_val;
    const vector<int> knapsack_val = { 360, 83, 59, 130, 431, 67, 230, 52, 93,
        125, 670, 892, 600, 38, 48, 147, 78, 256,
        63, 17, 120, 164, 432, 35, 92, 110, 22,
        42, 50, 323, 514, 28, 87, 73, 78, 15,
        26, 78, 210, 36, 85, 189, 274, 43, 33,
        10, 19, 389, 276, 312 };
    
    const vector<int> knapsack_weights = { 7, 0, 30, 22, 80, 94, 11, 81, 70,
        64, 59, 18, 0, 36, 3, 8, 15, 42,
        9, 0, 42, 47, 52, 32, 26, 48, 55,
        6, 29, 84, 2, 4, 18, 56, 7, 29,
        93, 44, 71, 3, 86, 66, 31, 65, 0,
        79, 20, 65, 52, 13 };
    //Knapsack recursive
    //cout << Knapsack_Recursive(KNAPSACK_SIZE - 1, KNAPSACK_CAPACITY, knapsack_val, knapsack_weights) << endl;
    
    //knapsack iterative/stack solution
    stack<tuple<int,int,int,int>> knapsack_stack; //1st capacity, 2nd array index, 3rd return value for max comparision,
    //4th stage for recursive call
    knapsack_stack.push(make_tuple(KNAPSACK_CAPACITY, KNAPSACK_SIZE-1, 0, 0));
    int ret_value = 0;
    int curridx, curr_capacity, curr_stage, curr_max_val;
    unordered_map<string,int> knapsack_memoize;
    while(!knapsack_stack.empty())
    {
        curr_max_val = get<2>(knapsack_stack.top());
        curridx = get<1>(knapsack_stack.top());
        curr_capacity = get<0>(knapsack_stack.top());
        curr_stage = get<3>(knapsack_stack.top());
        knapsack_stack.pop();
        
        string memoize_key = to_string(curridx) + ":" + to_string(curr_capacity);
        if(knapsack_memoize.find(memoize_key) != knapsack_memoize.end())
        {
            ret_value = knapsack_memoize[memoize_key];
            continue;
        }
        
        switch (curr_stage) {
            case 0:
            {
                if(0 == curr_capacity || curridx < 0)
                {
                    ret_value = 0;
                    knapsack_memoize.insert(make_pair(memoize_key, ret_value));
                }
                else if(knapsack_weights[curridx] > curr_capacity)
                {
                    //advance the stage to 3 and push it back to the stack (current snapshot needs to be processed
                    //after returning from recursive call
                    knapsack_stack.push(make_tuple(curr_capacity, curridx, 0, 3));
                    
                    //create a new stack element
                    knapsack_stack.push(make_tuple(curr_capacity, curridx-1, 0, 0));
                }
                else
                {
                    //advance the stage to 1 and push it back to the stack (current snapshot needs to be processed
                    //after returning from recursive call
                    knapsack_stack.push(make_tuple(curr_capacity, curridx, 0, 1));
                    
                    //create a new stack element
                    knapsack_stack.push(make_tuple(curr_capacity-knapsack_weights[curridx], curridx-1, 0, 0));
                }
                continue;
                break;
            }
            case 1:
            {
                //advance the stage to 2 and push it back to the stack (current snapshot needs to be processed
                //after returning from recursive call
                knapsack_stack.push(make_tuple(curr_capacity, curridx, ret_value+knapsack_val[curridx], 2));
                
                //create a new stack element
                knapsack_stack.push(make_tuple(curr_capacity, curridx-1, 0, 0));
                
                continue;
                break;
            }
            case 2:
            {
                ret_value = max(ret_value, curr_max_val);
                knapsack_memoize.insert(make_pair(memoize_key, ret_value));
                continue;
                break;
            }
            case 3:
            {
                knapsack_memoize.insert(make_pair(memoize_key, ret_value));
                continue;
                break;
            }
            default:
            {
                cerr << "unexpected stage " << curr_stage << endl;
                break;
            }
        }
    }
    */
    
    /*
    //stack/iterative solution
    stack<tuple<int, int, int, int>> knapsack_stack; //first - capacity, second - array index, third - return value for max comparision,
    //fourth - stage for recursive call
    knapsack_stack.push(make_tuple(KNAPSACK_CAPACITY, KNAPSACK_SIZE - 1, 0, 0));
    int ret_value = 0;
    int curridx, curr_capacity, curr_stage, curr_max_val;
    unordered_map<string,int> knapsack_memoize;
    while (!knapsack_stack.empty())
    {
        curr_max_val = get<2>(knapsack_stack.top());
        curridx = get<1>(knapsack_stack.top());
        curr_capacity = get<0>(knapsack_stack.top());
        curr_stage = get<3>(knapsack_stack.top());
        knapsack_stack.pop();
        
        string memoize_key = to_string(curridx) + ":" + to_string(curr_capacity);
        if (knapsack_memoize.find(memoize_key) != knapsack_memoize.end())
        {
            ret_value = knapsack_memoize[memoize_key];
            continue;
        }
        
        switch (curr_stage)
        {
            case 0:
            {
                if (0 == curr_capacity || curridx < 0)
                {
                    ret_value = 0;
                    knapsack_memoize.insert(make_pair(memoize_key, ret_value));
                }
                else if (knapsack_weights[curridx] > curr_capacity)
                {
                    //create a new stack element
                    knapsack_stack.push(make_tuple(curr_capacity, curridx - 1, 0, 0));
                }
                else
                {
                    //advance the stage to 1 and push it back in the stack (current snapshot needs to be processed after returning from recursive call
                    knapsack_stack.push(make_tuple(curr_capacity, curridx, 0, 1));
                    
                    //create a new stack element
                    knapsack_stack.push(make_tuple(curr_capacity - knapsack_weights[curridx], curridx - 1, 0, 0));
                }
                continue;
                break;
            }
            case 1:
            {
                //advance the stage to 2 and push it back in the stack (current snapshot needs to be processed after returning from recursive call
                knapsack_stack.push(make_tuple(curr_capacity, curridx, ret_value + knapsack_val[curridx], 2));
                
                //create a new stack element
                knapsack_stack.push(make_tuple(curr_capacity, curridx - 1, 0, 0));
                
                continue;
                break;
            }
            case 2:
            {
                ret_value = max(ret_value, curr_max_val);
                knapsack_memoize.insert(make_pair(memoize_key, ret_value));
                continue;
                break;
            }
            default:
            {
                cerr << "Unepxected stage" << endl;
                break;
            }
        }
    }
    */
    
    /*
    //Ugly numbers
    cout << "Isugly: " << IsUglyNum(21) << endl;
    
    //Print nth ugly number
    vector<int> UglyNums;
    const int UglyNumN = 40;
    int a2idx = 0, a3idx = 0, a5idx = 0;
    int a2val, a3val, a5val;
    
    cout << "1 ";
    UglyNums.push_back(1);
    for(auto i = 1; i < UglyNumN; ++i)
    {
        a2val = UglyNums[a2idx] * 2;
        a3val = UglyNums[a3idx] * 3;
        a5val = UglyNums[a5idx] * 5;
        const int min_val = min(a2val,min(a3val,a5val));
        UglyNums.push_back(min_val);
        if(min_val == a2val)
            ++a2idx;
        if(min_val == a3val)
            ++a3idx;
        if(min_val == a5val)
            ++a5idx;
        cout << min_val << " ";
    }
    cout << endl;
    */
    
    /*
    //Find max of sub-arrays of size k (k should be >= 2)
    const vector<int> input = {10,5,2,7,8,7,9,3,4,5,14,1,2,3,4,5,6,7,8,9};
    const int k = 8;
    int Max1 = max(input[0],input[1]), Max2 = min(input[0],input[1]);
    vector<int> KsubArrayOutput;
    for(auto i = 1; i < input.size(); ++i)
    {
        //Sliding window. Max1 and Max2 will hold the 1st and 2nd max elements
        //of the 'k sub-array'
        if(i < k)
        {
            Max1 = max(Max1,input[i]);
            if((Max2 < input[i]) && (input[i] < Max1))
                Max2 = input[i];
            if (i == k-1)
                KsubArrayOutput.push_back(Max1);
        }
        else
        {
            if(input[i-k] == Max1)
            {
                Max1 = max(Max2,input[i]);
                Max2 = min(Max2,input[i]);
            }
            else
            {
                int temp = Max1;
                Max1 = max(temp,input[i]);
                Max2 = min(temp,input[i]);
            }
            KsubArrayOutput.push_back(Max1);
        }
    }
    for(const auto &i : KsubArrayOutput)
        cout << i << " ";
    cout << endl;
    */
    
    /*
    //string to number conversion
    const string A = "1020";
    const int base = 16;
    cout << StrToNum(A,base) << endl;
     */
    
    /*
    //DFS for path traversal. True - wall, false - no wall
    //Min. steps to reach from start to end
    const int BoardN = 5, BoardM = 5;
    const pair<int,int> BoartStart = make_pair(4,4), BoardEnd = make_pair(0,0);
    queue<pair<int,int>> BoardBfsQ;
    const vector<vector<bool>> InBoard = {{false,true,true,true,false},
                                        {false,true,false,false,false},
                                        {false,true,false,true,false},
                                        {false,true,false,true,false},
                                        {false,false,false,false,false}};
    vector<vector<int>> BoardStepCnt (BoardM,vector<int>(BoardN,0));
    vector<vector<bool>> BoardVisited (BoardM,vector<bool>(BoardN,false));
    
    BoardBfsQ.push(BoardEnd);
    while(!BoardBfsQ.empty())
    {
        auto Curr = BoardBfsQ.front();
        BoardBfsQ.pop();
        BoardVisited[Curr.first][Curr.second] = true;
        
        if(Curr == BoartStart)
            break;
        
        //Add all 4 points
        auto Next = make_pair(Curr.first+1, Curr.second);
        if(IsBoardPntValid(Next, BoardM, BoardN, InBoard, BoardVisited))
        {
            BoardStepCnt[Next.first][Next.second] = BoardStepCnt[Curr.first][Curr.second] + 1;
            BoardBfsQ.push(Next);
        }
        Next = make_pair(Curr.first-1, Curr.second);
        if(IsBoardPntValid(Next, BoardM, BoardN, InBoard, BoardVisited))
        {
            BoardStepCnt[Next.first][Next.second] = BoardStepCnt[Curr.first][Curr.second] + 1;
            BoardBfsQ.push(Next);
        }
        Next = make_pair(Curr.first, Curr.second+1);
        if(IsBoardPntValid(Next, BoardM, BoardN, InBoard, BoardVisited))
        {
            BoardStepCnt[Next.first][Next.second] = BoardStepCnt[Curr.first][Curr.second] + 1;
            BoardBfsQ.push(Next);
        }
        Next = make_pair(Curr.first, Curr.second-1);
        if(IsBoardPntValid(Next, BoardM, BoardN, InBoard, BoardVisited))
        {
            BoardStepCnt[Next.first][Next.second] = BoardStepCnt[Curr.first][Curr.second] + 1;
            BoardBfsQ.push(Next);
        }
    }
    cout << BoardStepCnt[BoartStart.first][BoartStart.second] << endl;
    */
    
    /*
    //Check 2 strings which are max. 1 char diff - delete, insert or replace
    const string a = "abcd", b = "abc";
    const int LenDiff = static_cast<int>(a.length()) - static_cast<int>(b.length());
    if(abs(LenDiff) > 1)
    {
        cerr << "Error in len " << LenDiff << endl;
        return 0;
    }
    
    int i = 0, j = 0;
    bool IsFirstMismatch = false;
    
    while(i < a.length() && j < b.length())
    {
        if(a[i] != b[j])
        {
            if(!IsFirstMismatch)
            {
                IsFirstMismatch = true;
                if(1 == LenDiff)
                    ++i;
                else if(-1 == LenDiff)
                    ++j;
                else
                {
                    ++i;
                    ++j;
                }
                continue;
            }
            else
            {
                cerr << "1 char chk fail" << endl;
                return 0;
            }
        }
        else
        {
            ++i;
            ++j;
        }
    }
    cout << "1 char chk pass" << endl;
    */
    
    /*
    //Rotate NxN matrix to the right by 90deg
    const int MatRotSize = 5;
    vector<vector<int>> MatRot;
    for(auto i = 0; i < MatRotSize*MatRotSize; ++i)
    {
        if(0 == i % MatRotSize)
            MatRot.emplace_back(vector<int>());
        MatRot.back().push_back(i+1);
    }
    
    //Debug print before rotation
    for(auto i = 0; i < MatRotSize; ++i)
    {
        for(auto j = 0; j < MatRotSize; ++j)
            cout << MatRot[i][j] << " ";
        cout << endl;
    }
    
    int LoopCnt = 0, newN = MatRotSize;
    while(newN > 0)
    {
        int i = 0;
        while(i < newN-1)
        {
            //Swap 4 corner elements
            int Temp = MatRot[LoopCnt][LoopCnt+i];
            SwapMatElem(Temp, MatRot[LoopCnt+i][LoopCnt+newN-1]);
            SwapMatElem(Temp, MatRot[LoopCnt+newN-1][LoopCnt+newN-1-i]);
            SwapMatElem(Temp, MatRot[LoopCnt+newN-1-i][LoopCnt]);
            SwapMatElem(Temp, MatRot[LoopCnt][LoopCnt+i]);
            ++i;
        }
        newN -= 2;
        ++LoopCnt;
    }
        
    //Debug print after rotation
    cout << endl;
    for(auto i = 0; i < MatRotSize; ++i)
    {
        for(auto j = 0; j < MatRotSize; ++j)
            cout << MatRot[i][j] << " ";
        cout << endl;
    }
    */
    
    /*
    //Convert sorted array to BST
    const vector<int> SortedArray = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    BstHead = SortedArrayToBst(SortedArray, 0, SortedArray.size()-1);
    PrintBst();
     */
    
    /*
    //Identify number that repeats only once, all other integers in the input array repeat thrice
    const vector<int> NumRepeatThrice = {6,65535,3,6,3,6,3,19,19,19,13,13,13};
    vector<int> BitCnt (32,0); //assuming 32bit number
    for(const auto &i : NumRepeatThrice)
        ConvNumToBitsCnt(i,BitCnt);
    int RepeatOnceNum = 0;
    for(vector<int>::size_type i = 0; i < BitCnt.size(); ++i)
        if(BitCnt[i] % 3 != 0)
            RepeatOnceNum |= (1 << i);
    cout << "Num repeats only once " << RepeatOnceNum << endl;
    */
    
    /*
    //Print array in spiral pattern
    const int SpiralArrayM = 5, SpiralArrayN = 5;
    vector<vector<int>> SpiralInput;
    SpiralInput.resize(SpiralArrayM);
    int cnt = 0;
    for(auto i = 0; i < SpiralArrayM; ++i)
    {
        SpiralInput[i].resize(SpiralArrayN);
        for(auto j = 0; j < SpiralArrayN; ++j)
        {
            SpiralInput[i][j] = cnt++;
            cout << SpiralInput[i][j] << " ";
        }
        cout << endl;
    }
     
    int Top = 0, Bottom = SpiralArrayM-1, Left = 0, Right = SpiralArrayN-1;
    int Direction = 0; //0 - right, 1 - down, 2 - left, 3 - up
    while(Top <= Bottom && Left <= Right)
    {
        switch(Direction)
        {
            case 0:
                for(auto i = Left; i <= Right; ++i)
                    cout << SpiralInput[Top][i] << " ";
                Top++;
                break;
            case 1:
                for(auto i = Top; i <= Bottom; ++i)
                    cout << SpiralInput[i][Right] << " ";
                Right--;
                break;
            case 2:
                for(auto i = Right; i >= Left; --i)
                    cout << SpiralInput[Bottom][i] << " ";
                Bottom--;
                break;
            case 3:
                for(auto i = Bottom; i >= Top; --i)
                    cout << SpiralInput[i][Left] << " ";
                Left++;
                break;
        }
        ++Direction;
        Direction %= 4;
    }
    cout << endl;
    */
    
    /*
    //Find maximum value of contiguous subarray
    const vector<int> MaxSubarrayIn = {1,-3,4,-1,2,1};
    int MaxVal = 0, temp;
    for(auto i = 0; i < MaxSubarrayIn.size(); ++i)
        temp = FindArraySum(MaxSubarrayIn, i, MaxSubarrayIn.size()-1, MaxVal);
    cout << MaxVal << endl;
     */
    
    /*
    //Find maximum value of contiguous subarray
    const vector<int> MaxSubarrayIn = {1,-3,4,-1,2,1};
    int MaxHere, MaxSoFar;
    MaxHere = MaxSoFar = MaxSubarrayIn[0];
    for(int i = 1; i < MaxSubarrayIn.size(); ++i)
    {
        MaxHere = max(MaxHere + MaxSubarrayIn[i], MaxSubarrayIn[i]);
        MaxSoFar = max(MaxSoFar, MaxHere);
    }
    cout << MaxSoFar << endl;
    */
    
    /*
    //Num of paths from 1 corner to another - recursive
    //1 - path, 0 - no path. Can traverse only down or right
    //Find num paths from top left (0,0) to bottom right (M-1,N-1)
    #define SUM_PATHS_M    8
    #define SUM_PATHS_N    8
    const vector<vector<int>> SumPathsInput = { { 1,1,1,1,1,1,1,1 },
                                                { 1,1,0,1,1,1,0,1 },
                                                { 1,1,1,1,0,1,1,1 },
                                                { 0,1,0,1,1,0,1,1 },
                                                { 1,1,0,1,1,1,1,1 },
                                                { 1,1,1,0,0,1,0,1 },
                                                { 1,0,1,1,1,0,1,1 },
                                                { 1,1,1,1,1,1,1,1 } };
    cout << FindTotalPaths(SumPathsInput,SUM_PATHS_M,SUM_PATHS_N,0,0) << endl;
    */
    
    /*
    //Num of paths from 1 corner to another - iterative
    //1 - path, 0 - no path. Can traverse only down or right
    //Find num paths from top left (0,0) to bottom right (M-1,N-1)
    #define SUM_PATHS_M    8
    #define SUM_PATHS_N    8
    const vector<vector<int>> SumPathsInput = { { 1,1,1,1,1,1,1,1 },
                                                { 1,1,0,1,1,1,0,1 },
                                                { 1,1,1,1,0,1,1,1 },
                                                { 0,1,0,1,1,0,1,1 },
                                                { 1,1,0,1,1,1,1,1 },
                                                { 1,1,1,0,0,1,0,1 },
                                                { 1,0,1,1,1,0,1,1 },
                                                { 1,1,1,1,1,1,1,1 } };
    vector<vector<int>> SumPathsOutput (SUM_PATHS_M,vector<int>(SUM_PATHS_N,0));
    SumPathsOutput[SUM_PATHS_M-1][SUM_PATHS_N-1] = 1;
    for(int i = SUM_PATHS_M-1; i >= 0; --i)
        for(int j = SUM_PATHS_N-1; j>= 0; --j)
        {
            if((i == SUM_PATHS_M-1) && (j == SUM_PATHS_N-1))
                continue;
            if(!SumPathsInput[i][j])
                continue;
            int RightNgbr = 0;
            int BottomNgbr = 0;
            if((i+1 < SUM_PATHS_M) && (j < SUM_PATHS_N))
                BottomNgbr = SumPathsOutput[i+1][j];
            if((i < SUM_PATHS_M) && (j+1 < SUM_PATHS_N))
                RightNgbr = SumPathsOutput[i][j+1];
            SumPathsOutput[i][j] = BottomNgbr + RightNgbr;
        }
    cout << SumPathsOutput[0][0] << endl;
    */
    
    /*
    //Given array of inorder and preorder, build tree
    const vector<char> Inorder = {'d', 'b', 'e', 'a', 'f', 'c', 'g'};
    const vector<char> Preorder = {'a', 'b', 'd', 'e', 'c', 'f', 'g'};
    int StartIdx = 0;
    BstHead = BuildTreeFromInPreOrder(Inorder, Preorder, StartIdx, 0, Preorder.size()-1);
    PrintBst();
    */
    
    /*
    //LSD string sort (key indexed count)
    const vector<string> StringSortIn = {"4PGC938", "2IYE230", "3CIO720", "1ICK750", "1OHV845",
        "4JZY524", "1ICK750", "3CIO720", "*OHV845", "1OHV845", "2RLA629", "$RLA629", "3ATW723"};
    vector<string> StringSortOut = StringSortIn;
    vector<string> TempSort (StringSortOut.size()," ");
    for(int i = StringSortOut[0].size()-1; i >= 0; --i) //for every character in string (all string should have equal length)
    {
        const int MaxChar = 128; //ASCII
        vector<int> CharCnt (MaxChar+1,0); //Have char cnt size + 1 of total num of different characters
        
        //Cnt characters
        for(int j = 0; j < StringSortOut.size(); ++j) //for all string inputs
            CharCnt[static_cast<int>(StringSortOut[j][i]) + 1]++;
        
        //Convert cnt to indices
        for(int j = 0; j < MaxChar; ++j) //for char cnt array
            CharCnt[j+1] += CharCnt[j];
        
        //Build temp array based on char cnt index (stable sort)
        for(int j = 0; j < TempSort.size(); ++j) //for all string inputs
            TempSort[ CharCnt[static_cast<int>(StringSortOut[j][i])]++ ] = StringSortOut[j];
        
        //Copy the array back
        for(int j = 0; j < TempSort.size(); ++j) //for all string inputs
            StringSortOut[j] = TempSort[j];
    }
    
    //Print unsorted and sorted output
    for(int j = 0; j < StringSortOut.size(); ++j) //for all string inputs
        cout << StringSortIn[j] << "  " << StringSortOut[j] << endl;
    */
    
    /*
    //check if array of integers can be split into 2 subarrays with equal sum
    const vector<int> SplitSumInput {-5,5,-5,-5,-5};
    cout << IsSplitArraySumPossible(SplitSumInput, SplitSumInput.size()-1, true, 0, 0) << endl;
    */
    
    /*
    //find exponential i.e. x^y power(x,y)
    cout << CustomPower(2,10) << endl;
    */
    
    /*
    //DCP69 - largest product that can be made by multiplying any three integers.
    //Below solution doesnt work for some combination of -ve numbers
    const vector<int> ProdArray {13,10,5,-200,0,-200};
    cout << MaxProduct(ProdArray, ProdArray.size()-1, 3) << endl;
    */
    
    /*
    //DCP70 - nth perfect number - number that sums to 10
    //every perfect number is a multiple of 9 + 1. 1st number is 19.
    //so start with 19 and increment by 9
    int PerfectNumCnt = 0;
    const int N = 144;
    int i = 19;
    while(1)
    {
        int temp = i;
        int sum = 0;
        while(temp != 0)
        {
            sum += (temp % 10);
            temp /= 10;
        }
        if(10 == sum)
        {
            ++PerfectNumCnt;
            cout << i << " " << PerfectNumCnt << endl;
            if(PerfectNumCnt == N)
            {
                cout << N << "th perfect number is " << i << endl;
                break;
            }
        }
        i += 9;
    }
     */
    
    /*
    //Number of bishops attack
    const vector<pair<int,int>> BishopsAttackIn {{make_pair(0,0)}, {make_pair(2,2)}, {make_pair(4,0)}, {make_pair(4,4)}, {make_pair(0,3)},
        {make_pair(1,4)}};
    const int BishopsChessSize = 5;
    vector<vector<int>> BishopsTemp (BishopsChessSize, vector<int>(BishopsChessSize,0));
    int AttackCnt = 0;
    for(int i = 0; i < BishopsAttackIn.size(); ++i)
    {
        AttackCnt += BishopsTemp[BishopsAttackIn[i].first][BishopsAttackIn[i].second];
        LogBishopPath(BishopsTemp, BishopsAttackIn[i].first, BishopsAttackIn[i].second, BishopsChessSize);
    }
    cout << AttackCnt << endl;
    */
    
    /*
    //Overlapping intervals (DCP77)
    vector<pair<int,int>> Intervals {make_pair(1, 3), make_pair(5, 8), make_pair(4, 10), make_pair(20, 25), make_pair(1,10)};
    vector<pair<int,int>> OutputIntervals;
    //Sort the input based on start times (first element in pair)
    sort(Intervals.begin(), Intervals.end(), PairSortComp);
    pair<int,int> MinMax = make_pair(INT_MAX, INT_MIN);
    
    for(int i = 0; i < Intervals.size(); ++i)
    {
        if(Intervals[i].second > MinMax.second)
        {
            MinMax.second = Intervals[i].second;
            OutputIntervals.push_back(Intervals[i]);
        }
    }
    */
    
    /*
    //DCP #84 - return number of islands of 1
    #define ISLANDS_M    6
    #define ISLANDS_N    5
    const vector<vector<int>> IslandsInput = { { 1,0,0,0,0 },
                                               { 0,0,1,1,0 },
                                               { 0,1,1,0,1 },
                                               { 0,0,0,0,0 },
                                               { 1,1,0,0,1 },
                                               { 1,1,0,0,1 } };
    int IslandsCnt = 0;
    vector<vector<bool>> Visited {ISLANDS_M, vector<bool>(ISLANDS_N, false)};
    for(int row = 0; row < ISLANDS_M; ++row)
    {
        for(int col = 0; col < ISLANDS_N; ++col)
        {
            if(!Visited[row][col])
            {
                Visited[row][col] = true;
                if(IslandsInput[row][col])
                {
                    IslandsBfs(row, col, Visited, IslandsInput);
                    ++IslandsCnt;
                }
            }
        }
    }
    */
    
    /*
    //Find the length of Longest Increasing Subsequence - LIS (DCP #75)
    //const vector<int> InputVec = {10,8,3,9,2,11,4,7,12,13,14};
    const vector<int> InputVec = {6,4,5,3,2,1,7,8,9,10,11};
    vector<int> LisArray (static_cast<int>(InputVec.size()), 1);
    int MaxLis = LisArray[0];
    for(int i = 1; i < InputVec.size(); ++i)
    {
        for(int j = 0; j < i; ++j)
        {
            if(InputVec[i] > InputVec[j])
                LisArray[i] = max(LisArray[i], LisArray[j]+1);
        }
        if(MaxLis < LisArray[i])
            MaxLis = LisArray[i];
    }
    */
    
    /*
    //largest sum of non-adjacent numbers. DCP #9
    const vector<int> InputVec = {2,-400,6,200,5};
    unordered_map<int,int> Memo;
    cout << LargestNonAdjSum(InputVec, InputVec.size()-1, Memo) << endl;
    */
    
    /*
    //Running median
    priority_queue<int, vector<int>, greater<int>> RunMedMinHeap;
    priority_queue<int> RunMedMaxHeap;
    const vector<int> InputVec = {6,4,5,3,2,1,7,8,9,10,11};
    for(const auto &i : InputVec)
    {
        //Split elements into min. heap and max. heap
        //all elements <= median should be in max heap.
        //all elements >= median should in min heap
        AddElemMedian(i, RunMedMinHeap, RunMedMaxHeap);
        BalanceMedian(RunMedMinHeap, RunMedMaxHeap);
        cout << GetMedian(RunMedMinHeap, RunMedMaxHeap) << endl;
    }
    */
    
    /*
    //All permutations of a string
    const string input_str {"abcde"};
    unordered_set<string> StrPermOutputMap;
    for(string::size_type i = 0; i < input_str.size(); ++i)
    {
        vector<string> temp_vec;
        temp_vec.push_back(input_str.substr(i,1));
        StrPermOutputMap.insert(input_str.substr(i,1));
        for(string::size_type j = i+1; j < input_str.size(); ++j)
            StringPerm(temp_vec, j, input_str, StrPermOutputMap);
    }
    
    map<int, int> SizeSet;
    for(string::size_type i = 0; i < input_str.size(); ++i)
        SizeSet.insert(make_pair(i+1, 0));
    for(const auto &i : StrPermOutputMap)
    {
        cout << i << endl;
        ++SizeSet[i.size()];
    }
    for(const auto &i : SizeSet)
        cout << i.second << " ";
    cout << endl;
    */
    
    return 0;
}

void StringPerm(const vector<string> &Input, const string::size_type idx, const string &input_str, unordered_set<string> &StrPermOutputMap)
{
    if(idx >= input_str.size())
        return;
    vector<string> temp_vec;
    string temp;
    for(string::size_type i = 0; i < Input.size(); ++i)
    {
        const int StrLen = Input[i].size();
        temp = input_str[idx] + Input[i];
        
        if(StrPermOutputMap.find(temp) != StrPermOutputMap.end())
            continue; //permutation already identified
        
        temp_vec.push_back(temp);
        StrPermOutputMap.insert(temp);
        for(string::size_type k = 0; k < Input[i].size(); ++k)
        {
            temp = Input[i].substr(0,k+1) + input_str.substr(idx,1) + Input[i].substr(k+1, StrLen-k-1);
            if(StrPermOutputMap.find(temp) != StrPermOutputMap.end())
                continue; //permutation already identified
            
            temp_vec.push_back(temp);
            StrPermOutputMap.insert(temp);
        }
        for(string::size_type k = idx+1; k < input_str.size(); ++k)
            StringPerm(temp_vec, k, input_str, StrPermOutputMap);
    }
}

void AddElemMedian(const int i, priority_queue<int, vector<int>, greater<int>> &MinHeap, priority_queue<int> &MaxHeap)
{
    if(MaxHeap.empty())
        MaxHeap.push(i);
    else
    {
        if(i < MaxHeap.top())
            MaxHeap.push(i);
        else
            MinHeap.push(i);
    }
}

void BalanceMedian(priority_queue<int, vector<int>, greater<int>> &MinHeap, priority_queue<int> &MaxHeap)
{
    if((MinHeap.size() - MaxHeap.size()) == 2)
    {
        int temp = MinHeap.top();
        MinHeap.pop();
        MaxHeap.push(temp);
    }
    else if ((MaxHeap.size() - MinHeap.size()) == 2)
    {
        int temp = MaxHeap.top();
        MaxHeap.pop();
        MinHeap.push(temp);
    }
    else if (abs(static_cast<int>(MaxHeap.size()) - static_cast<int>(MinHeap.size())) > 2)
    {
        cout << "Error!! Something went wrong. MinHeap size " << MinHeap.size() << " MaxHeap size " << MaxHeap.size() << endl;
    }
}

float GetMedian(const priority_queue<int, vector<int>, greater<int>> &MinHeap, const priority_queue<int> &MaxHeap)
{
    if(MinHeap.size() == MaxHeap.size())
        return( static_cast<float>(MinHeap.top()+MaxHeap.top())/2.0 );
    return( MinHeap.size() > MaxHeap.size() ? MinHeap.top() : MaxHeap.top() );
}

int LargestNonAdjSum(const vector<int> &InputVec, const int idx, unordered_map<int,int> &Memo)
{
    if(idx < 0)
        return 0;
    if(0 == idx)
        return InputVec[0];
    if(Memo.find(idx) != Memo.end())
        return Memo[idx];
    
    int ret = max(LargestNonAdjSum(InputVec, idx-2, Memo) + InputVec[idx],
                  LargestNonAdjSum(InputVec, idx-1, Memo));
    Memo.insert(make_pair(idx,ret));
                  
    return(ret);
}

bool IsIslandValid(const pair<int,int> NextPnt, const int M, const int N,
                   const vector<vector<int>> &IslandsInput,
                   const vector<vector<bool>> &Visited)
{
    if(NextPnt.first < 0 || NextPnt.first >= M || NextPnt.second < 0 || NextPnt.second >= N)
        return false;
    if(!IslandsInput[NextPnt.first][NextPnt.second] || Visited[NextPnt.first][NextPnt.second])
        return false;
    return true;
}

void IslandsBfs(const int row, const int col, vector<vector<bool>> &Visited, const vector<vector<int>> &IslandsInput)
{
    queue<pair<int,int>> NgbrQ;
    NgbrQ.push(make_pair(row, col));
    while(!NgbrQ.empty())
    {
        pair<int,int> Curr = NgbrQ.front(), Next;
        NgbrQ.pop();
        Visited[Curr.first][Curr.second] = true;
        
        Next = make_pair(Curr.first+1, Curr.second);
        if(IsIslandValid(Next, Visited.size(), Visited[0].size(), IslandsInput, Visited))
            NgbrQ.push(Next);
        
        Next = make_pair(Curr.first, Curr.second+1);
        if(IsIslandValid(Next, Visited.size(), Visited[0].size(), IslandsInput, Visited))
            NgbrQ.push(Next);
        
        Next = make_pair(Curr.first-1, Curr.second);
        if(IsIslandValid(Next, Visited.size(), Visited[0].size(), IslandsInput, Visited))
            NgbrQ.push(Next);
        
        Next = make_pair(Curr.first, Curr.second-1);
        if(IsIslandValid(Next, Visited.size(), Visited[0].size(), IslandsInput, Visited))
            NgbrQ.push(Next);
    }
}
    
bool PairSortComp(pair<int,int> i, pair<int,int> j)
{
    if(i.first == j.first)
        return (j.second < i.second);
    return (i.first < j.first);
}

bool CHECK_BISHOP_VALID_CHESS_CELL(const int R, const int C, const int N)
{
    if(R<0 || R>=N)
        return false;
    if(C<0 || C>=N)
        return false;
    return true;
}
void LogBishopPath(vector<vector<int>> &Temp, const int R, const int C, const int N)
{
    for(int i = 0; i < N; ++i)
    {
        if(CHECK_BISHOP_VALID_CHESS_CELL(R+i,C+i,N))
            Temp[R+i][C+i] = 1;
        if(CHECK_BISHOP_VALID_CHESS_CELL(R-i,C+i,N))
            Temp[R-i][C+i] = 1;
        if(CHECK_BISHOP_VALID_CHESS_CELL(R+i,C-i,N))
            Temp[R+i][C-i] = 1;
        if(CHECK_BISHOP_VALID_CHESS_CELL(R-i,C-i,N))
            Temp[R-i][C-i] = 1;
    }
}

int MaxProduct(const vector<int> &A, const int idx, const int cnt)
{
    static int loopcnt = 0;
    ++loopcnt;
    int ret;
    
    cout << "entry " << idx << " " << cnt << " " << endl;
    
    if(( idx < 0 ) || (cnt <= 0 ))
        ret = 1;
    else if(idx < cnt)
        ret = MaxProduct(A, idx-1, cnt-1) * A[idx];
    else
        ret = max( MaxProduct(A, idx-1, cnt), MaxProduct(A, idx-1, cnt-1)*A[idx] );
    
    cout << "ret " << ret << " " << idx << " " << cnt << " " << endl;
    
    return ret;
}

long int CustomPower(const int x, const int y)
{
    static int LoopCtr = 0;
    ++LoopCtr;
    
    static unordered_map<int, long int> PowerMap;
    if (y <= 0)
        return 1;
    if (1 == y)
        return x;
    if(PowerMap.find(y) != PowerMap.end())
        return PowerMap[y];
    
    long int ret;
    if(y % 2)
        ret = CustomPower(x,y/2) * CustomPower(x,y/2+1); //odd y
    else
        ret = CustomPower(x,y/2) * CustomPower(x,y/2); //even y
    
    PowerMap.insert(make_pair(y,ret));
    
    cout << ret << " " << x << " " << y << " " << LoopCtr << endl;
    
    return ret;
}

bool IsSplitArraySumPossible(const vector<int> &A, const int idx, const bool IsB1, int B1sum, int B2sum)
{
    static int LoopCtr = 0;
    ++LoopCtr;
    
    if(IsB1)
        B1sum += A[idx];
    else
        B2sum += A[idx];
 
    cout << idx << " " << IsB1 << " " << B1sum << " " << B2sum << " " << LoopCtr << endl;
    
    if(idx <= 0)
        return (B1sum == B2sum);
    
    return (IsSplitArraySumPossible(A, idx-1, true, B1sum, B2sum) ||
            IsSplitArraySumPossible(A, idx-1, false, B1sum, B2sum));
}

shared_ptr<BstType> BuildTreeFromInPreOrder
(const vector<char> &Inorder, const vector<char> &Preorder, int &PreorderIdx, const int lo, const int hi)
{
    cout << "BuildTreeFromInPreOrder: lo " << lo << " hi " << hi << " PreIdx " << PreorderIdx << endl;
    const auto N = Preorder.size();
    if(lo > hi || PreorderIdx >= N || lo < 0 || hi >= N)
        return nullptr;
    
    int InorderIdx = FindInorderIdx(Inorder, Preorder[PreorderIdx], lo, hi);
    if(-1 == InorderIdx)
        return nullptr;
    
    shared_ptr<BstType> NewNode = make_shared<BstType>();
    NewNode->charkey = Preorder[PreorderIdx];
    PreorderIdx++;
    NewNode->left = BuildTreeFromInPreOrder(Inorder, Preorder, PreorderIdx, lo, InorderIdx-1);
    NewNode->right = BuildTreeFromInPreOrder(Inorder, Preorder, PreorderIdx, InorderIdx+1, hi);
    NewNode->size = GetSize(NewNode->left) + GetSize(NewNode->right) + 1;
    
    return NewNode;
}

int FindInorderIdx(const vector<char> &Inorder, const char val, const int lo, const int hi)
{
    cout << "FindInorderIdx: lo " << lo << " hi " << hi << " val " << val << endl;
    
    if(lo > hi || lo < 0 || hi >= Inorder.size())
        return -1;
    
    for(auto i = lo; i <= hi; ++i)
        if(Inorder[i] == val)
            return i;
    
    return -1;
}

int FindTotalPaths(const vector<vector<int>> &SumPathsInput, const int M, const int N, const int row, const int col)
{
    static vector<vector<int>> Memo(M,vector<int>(N,-1));
    if((row >= M) || (col >= N))
        return 0;
    if(0 == SumPathsInput[row][col])
        return 0;
    if(-1 != Memo[row][col])
        return Memo[row][col];
    if((row == M-1) && (col == N-1))
        return 1;
    Memo[row][col] = FindTotalPaths(SumPathsInput,M,N,row+1,col) +
        FindTotalPaths(SumPathsInput,M,N,row,col+1);
    return Memo[row][col];
}

int FindArraySum(const vector<int> &A, const int Start, const int End, int &MaxVal)
{
    static unordered_map<string,int> ArraySumMap;
    static int LoopCnt = 0;
    if(Start > End)
        return 0;
    
    string IdxInput = to_string(Start) + ":" + to_string(End);
    if(ArraySumMap.find(IdxInput) != ArraySumMap.end())
        return ArraySumMap[IdxInput];
    
    LoopCnt++;
    cout << "LoopCnt " << LoopCnt << " Start " << Start << " End " << End << endl;
    
    int ret = A[End] + FindArraySum(A,Start,End-1,MaxVal);
    if(ret > MaxVal)
        MaxVal = ret;
    ArraySumMap.insert(make_pair(IdxInput,ret));
    return ret;
}

void ConvNumToBitsCnt(const int num, vector<int> &BitCnt)
{
    int i = 0;
    int mask = 1 << i;
    while (mask <= num)
    {
        if(mask & num)
            ++BitCnt[i];
        ++i;
        mask = 1 << i;
    }
}

shared_ptr<BstType> SortedArrayToBst(const vector<int> A, const int lo, const int hi)
{
    cout << "SortedArrayToBst: lo " << lo << " hi " << hi << endl;
    if(lo > hi)
        return nullptr;
    int mid = (lo+hi)/2;
    shared_ptr<BstType> NewNode = make_shared<BstType>();
    NewNode->key = A[mid];
    NewNode->left = SortedArrayToBst(A, lo, mid-1);
    NewNode->right = SortedArrayToBst(A, mid+1, hi);
    NewNode->size = GetSize(NewNode->left) + GetSize(NewNode->right) + 1;
    return NewNode;
}

void SwapMatElem(int &Val1, int &Val2)
{
    int Temp = Val2;
    Val2 = Val1;
    Val1 = Temp;
}

bool IsBoardPntValid(const pair<int,int> NextPnt, const int M, const int N,
                     const vector<vector<bool>> &InBoard,
                     const vector<vector<bool>> &BoardVisited)
{
    if(NextPnt.first < 0 || NextPnt.first >= M || NextPnt.second < 0 || NextPnt.second >= N)
        return false;
    if(InBoard[NextPnt.first][NextPnt.second] || BoardVisited[NextPnt.first][NextPnt.second])
        return false;
    return true;
}

int StrToNum(const string &In, const int base)
{
    int RetVal = 0;
    int Mantissa = 0;
    auto StrLen = In.length();
    for(string::size_type i = 0; i < StrLen; ++i)
    {
        Mantissa = stoi(In.substr(i,1));
        if(Mantissa >= base)
        {
            cerr << "Invalid num " << Mantissa << " Base " << base << endl;
            return -99;
        }
        RetVal += (Mantissa * static_cast<int>(pow(base,StrLen-i-1)));
    }
    return RetVal;
}

const vector<int> UglyNumFactors = {2,3,5};
bool IsUglyNum(const int n)
{
    bool ret = false;
    if(n <= 0)
        return false;
    if(1 == n)
        return true;
    for(const auto &i : UglyNumFactors)
    {
        if(0 == n % i)
            ret = IsUglyNum(n/i);
        if(ret)
            return true;
    }
    return false;
}

unordered_map<string,int> knapsack_memoize_rec;
int Knapsack_Recursive(const int curridx, const int curr_capacity, const vector<int> &knapsack_val, const vector<int> &knapsack_weights)
{
    static int loopcnt = 0;
    cout << "Rec fn entry: " << curr_capacity << " " << curridx << " " << ++loopcnt << endl;
    
    int ret;
    string memoize_key = to_string(curridx) + ":" + to_string(curr_capacity);
    if (knapsack_memoize_rec.find(memoize_key) == knapsack_memoize_rec.end())
    {
        if (0 == curr_capacity || curridx < 0)
        {
            cout << "ret 0." << endl;;
            ret = 0;
        }
        else if (knapsack_weights[curridx] > curr_capacity)
        {
            cout << "wight curr idx > capacity." << endl;
            ret = Knapsack_Recursive(curridx - 1, curr_capacity, knapsack_val, knapsack_weights);
        }
        else
        {
            cout << "Yes or no case." << endl;
            int temp1, temp2;
            temp1 = Knapsack_Recursive(curridx - 1, curr_capacity - knapsack_weights[curridx], knapsack_val, knapsack_weights) + knapsack_val[curridx];
            temp2 = Knapsack_Recursive(curridx - 1, curr_capacity, knapsack_val, knapsack_weights);
            ret = max(temp1, temp2);
        }
        knapsack_memoize_rec.insert(make_pair(memoize_key,ret));
        cout << "              Return value: " << ret << " " << curr_capacity << " " << curridx << endl;
        return ret;
    }
    else
    {
        cout << "Return memoize: " << knapsack_memoize_rec[memoize_key] << endl;
        return knapsack_memoize_rec[memoize_key];
    }
}

bool NqueensProblemSol2(const int Col)
{
    NqueensLoopCtr++;
    
    if (Col >= CHESS_BOARD_SIZE)
        return true;
    int Row = 0;
    while (Row < CHESS_BOARD_SIZE)
    {
        NqueensLocVec.push_back(make_pair(Row, Col));
        if (IsQueenPosValid(Row, Col) && NqueensProblemSol2(Col + 1))
            return true;
        NqueensLocVec.pop_back();
        ++Row;
    }
    return false;
}

bool IsQueenPosValid(const int Row, const int Col)
{
    //Compare only until 1 minus the last element
    for (int k = 0; (NqueensLocVec.size() > 1) && (k < NqueensLocVec.size()-1); ++k)
    {
        //Compare only rows and diagonal elements. Column need not be compared
        //as we are filling-up the chess board from left to right
        if ((Row == NqueensLocVec[k].first) ||
            (abs(Row - NqueensLocVec[k].first) == abs(Col - NqueensLocVec[k].second)))
            return false;
    }
    return true;
}

void NqueensProblem(void)
{
    NqueensLoopCtr++;
    
    for (int i = 0; i < CHESS_BOARD_SIZE; ++i)
        for (int j = 0; j < CHESS_BOARD_SIZE; ++j)
            if (NextQueenLocTable[i][j])
            {
                NqueensLocVec.push_back(make_pair(i, j));
                //Update possible next queen locations array
                UpdateNextQueenLocTable();
                NqueensProblem();
            }
    
    //At this point, there are no more possible squares for Queen. Check total size
    if (NqueensLocVec.size() != CHESS_BOARD_SIZE)
    {
        //Remove the previous queen position and iterate again
        NqueensLocVec.pop_back();
        UpdateNextQueenLocTable();
    }
}

void UpdateNextQueenLocTable(void)
{
    bool CollisionFound = false;
    for (int i = 0; i < CHESS_BOARD_SIZE; ++i)
        for (int j = 0; j < CHESS_BOARD_SIZE; ++j)
        {
            CollisionFound = false;
            for(auto k = NqueensLocVec.begin(); k != NqueensLocVec.end(); ++k)
            {
                if ((i == k->first) || (j == k->second) || (abs(i - k->first) == abs(j - k->second)))
                {
                    CollisionFound = true;
                    NextQueenLocTable[i][j] = false;
                    break;
                }
            }
            if (!CollisionFound)
                NextQueenLocTable[i][j] = true;
        }
}

bool IsValidChessSquare(const int R, const int C)
{
    return ((R >= 0) && (R < CHESS_BOARD_SIZE) &&
            (C >= 0) && (C < CHESS_BOARD_SIZE));
}

void KnightsTourAlgm(const int R, const int C)
{
    int CurrBoardSq = CHESS_RC_TO_BOARD(R, C);
    SqVisited.set(CurrBoardSq);
    if (SqVisited.all())
        LastSqVisited = CurrBoardSq;
    PreTraverse.push_back(make_pair(R, C));
    for (const auto &i : KnightDirections)
    {
        pair<int, int> NextSq = make_pair(R + i.first, C + i.second);
        int NextBoardSq = CHESS_RC_TO_BOARD(NextSq.first, NextSq.second);
        if (IsValidChessSquare(NextSq.first, NextSq.second))
        {
            if (!SqVisited.test(NextBoardSq))
            {
                SqParent[NextBoardSq] = CurrBoardSq;
                KnightsTourAlgm(NextSq.first, NextSq.second);
            }
        }
    }
    PostTraverse.push_back(make_pair(R, C));
}

void InsertKey(const int key, const string value)
{
    BstHead = InsertKey(BstHead, key, value);
    cout << "Insert complete: Key: " << key << " Value: " << value << endl;
}

shared_ptr<BstType> InsertKey(shared_ptr<BstType> Node, const int key, const string value)
{
    //Key is not present. Create a new node, and update size
    if (nullptr == Node)
    {
        shared_ptr<BstType> NewNode = make_shared<BstType>();
        NewNode->key = key;
        NewNode->value = value;
        return NewNode;
    }
    
    if (key == Node->key)
        //Key is already present. Just replace value with new one and return;
        Node->value = value;
    else if (key > Node->key)
        Node->right = InsertKey(Node->right, key, value);
    else
        Node->left = InsertKey(Node->left, key, value);
    
    Node->size = GetSize(Node->left) + GetSize(Node->right) + 1;
    return Node;
}

void DeleteKey(const int key)
{
    BstHead = DeleteKey(BstHead, key);
}

shared_ptr<BstType> DeleteKey(shared_ptr<BstType> Node, const int key)
{
    if (nullptr == Node)
    {
        cout << "Delete: Key: " << key << " not found" << endl;
        return Node;
    }
    
    if (key > Node->key)
        Node->right = DeleteKey(Node->right, key);
    else if (key < Node->key)
        Node->left = DeleteKey(Node->left, key);
    else
    {
        cout << "Delete: Key: " << key << endl;
        
        //Key found
        if (nullptr == Node->left)
        {
            shared_ptr<BstType> ret = Node->right;
            Node.reset();
            return ret;
        }
        else if (nullptr == Node->right)
        {
            shared_ptr<BstType> ret = Node->left;
            Node.reset();
            return ret;
        }
        else
        {
            //Since key to be deleted has 2 childs, Replace the node to be deleted with its min. node of right tree
            shared_ptr<BstType> MinNode = GetMin(Node->right);
            Node->key = MinNode->key;
            Node->value = MinNode->value;
            Node->right = DeleteMin(Node->right); //Delete the min. node of right tree
        }
    }
    
    Node->size = GetSize(Node->left) + GetSize(Node->right) + 1;
    return Node;
}

string GetValue(const int key)
{
    shared_ptr<BstType> temp = BstHead;
    while (nullptr!=temp)
    {
        if (key == temp->key)
            return temp->value;
        else if (key > temp->key)
            temp = temp->right;
        else
            temp = temp->left;
    }
    return "Not Found!!!";
}

void DeleteBst(void)
{
    DeleteBst(BstHead);
    BstHead = nullptr;
    cout << "Delete BST complete" << endl;
}

void DeleteBst(shared_ptr<BstType> Node)
{
    if (nullptr == Node)
        return;
    DeleteBst(Node->left);
    DeleteBst(Node->right);
    Node.reset();
}

void PrintBst(void)
{
    PrintBst(BstHead);
    cout << "Print BST complete" << endl;
}

void PrintBst(const shared_ptr<BstType> Node)
{
    if (nullptr == Node)
        return;
    PrintBst(Node->left);
    cout << "Print: Key: " << Node->key << " Value: " << Node->value << " Char: " <<  Node->charkey << " Size: " << Node->size << endl;
    PrintBst(Node->right);
}

int GetSize(const shared_ptr<BstType> Node)
{
    if (nullptr == Node)
        return 0;
    return Node->size;
}

int GetRank(const int key)
{
    return GetRank(BstHead, key);
}

int GetRank(const shared_ptr<BstType> Node, const int key)
{
    if (nullptr == Node)
        return 0;
    
    if (key == Node->key)
        return GetSize(Node->left);
    else if (key < Node->key)
        return GetRank(Node->left, key);
    else
    {
        return (GetRank(Node->right, key) + GetSize(Node->left)+1);
    }
}

int GetMin(void)
{
    if (nullptr != BstHead)
        return GetMin(BstHead)->key;
    else
        return -9999;
}

shared_ptr<BstType> GetMin(const shared_ptr<BstType> Node)
{
    if (nullptr == Node->left)
        return Node;
    return GetMin(Node->left);
}

void DeleteMin(void)
{
    if(nullptr != BstHead)
        BstHead = DeleteMin(BstHead);
}

shared_ptr<BstType> DeleteMin(shared_ptr<BstType> Node)
{
    if (nullptr == Node->left)
    {
        shared_ptr<BstType> ret = Node->right;
        Node.reset();
        return ret;
    }
    Node->left = DeleteMin(Node->left);
    Node->size = GetSize(Node->left) + GetSize(Node->right) + 1;
    return Node;
}

int GetMax(void)
{
    return 0;
}

int NumWaysNsum(const int Sum, const vector<int>::size_type VecIdx, unordered_map<string,int> &ResultMap, const vector<int> &InputVec)
{
    static int FuncCnt = 0;
    
    ++FuncCnt;
    cout << "Total: " << Sum << " Inputs: ";
    for(auto i = VecIdx; i < InputVec.size(); --i)
        cout << InputVec[i] << " ";
    cout << "FuncCnt: " << FuncCnt << endl;
    
    if(0 == Sum)
        return 1;
    else if (Sum < 0)
        return 0;
    else if(VecIdx >= InputVec.size())
        return 0;
    
    int Ret = 0;
    string key = to_string(Sum) + ':' + to_string(VecIdx);
    auto MapItr = ResultMap.find(key);
    if(MapItr != ResultMap.end())
        return MapItr->second;
    
    if(Sum < InputVec[VecIdx])
        Ret = NumWaysNsum(Sum, VecIdx-1, ResultMap, InputVec);
    else
        Ret = NumWaysNsum(Sum-InputVec[VecIdx], VecIdx-1, ResultMap, InputVec) + NumWaysNsum(Sum, VecIdx-1, ResultMap, InputVec);
    ResultMap[key] = Ret; //Create new element in Map and store the result
    
    cout << "Ret: Total: " << Sum << " Inputs: ";
    for(auto i = VecIdx; i < InputVec.size(); --i)
        cout << InputVec[i] << " ";
    cout << "NumWays: " << Ret << endl;
    
    return Ret;
}

int NumWaysStair(const int Total, unordered_map<int,int> &ResultMap, const vector<int> &Steps)
{
    static int FuncCnt = 0;
    
    if(0 == Total)
        return 1;
    
    auto MapItr = ResultMap.find(Total);
    if(MapItr != ResultMap.end())
        return MapItr->second;
    
    ++FuncCnt;
    cout << "Total: " << Total << " FuncCnt: " << FuncCnt << endl;
    
    int Ret = 0;
    for(auto &i:Steps)
        if(Total-i>=0)
            Ret += NumWaysStair(Total-i, ResultMap, Steps);
    
    ResultMap[Total] = Ret; //Create new element in Map and store the result
    
    cout << "Ret: Total: " << Total << " NumWays: " << Ret << endl;
    
    return Ret;
}

/*
 void SetRoot(const int R, const int C, const int RootVal)
 {
 if((R<0) || (R>=M)||
 (C<0) || (C>=N))
 return;
 if((RootVal == ComponentRoot[R][C]) || (!ComponentInputs[R][C]))
 return;
 
 //Set RootVal for current node and for all its 4 ngbrs
 ComponentRoot[R][C] = RootVal;
 for(auto &k : MatrixNgbrs)
 SetRoot(R+k.first, C+k.second, RootVal);
 }
 */

void InsertHeapTree(vector<int> &heaptree, int val)
{
    //Add new value to the end of the tree, and then swim up to maintain heap
    heaptree.push_back(val);
    HeapSwimUp(heaptree, heaptree.size()-1);
}

int GetHeapTree(vector<int> &heaptree)
{
    //return 0th index which will be the max. After returning the 0th index,
    //copy the last element into 0th index, reduce the tree size and then
    //swim down to maintain heap
    if(!heaptree.size())
        return -9999;
    
    int ret = heaptree[0];
    SwapIvec(heaptree,0,heaptree.size()-1);
    heaptree.pop_back();
    if(heaptree.size())
        HeapSwimDown(heaptree,0,heaptree.size());
    
    return ret;
}

void HeapSwimUp(vector<int> &heaptree, vector<int>::size_type k)
{
    //Parent of k will be (k-1)/2
    //Keep swapping until parent > child
    while((k>0) && (k<heaptree.size()))
    {
        if(heaptree[(k-1)/2] < heaptree[k])
        {
            SwapIvec(heaptree,k,(k-1)/2);
            k = (k-1)/2;
        }
        else
            break;
    }
}

void HeapSwimDown(vector<int> &heaptree, vector<int>::size_type k, vector<int>::size_type vecsize)
{
    //Children of k will be (2k+1) and (2K+2)
    //Keep swapping until parent > max(child1,child2)
    while(2*k+1 < vecsize)
    {
        auto c1 = 2*k+1, max = c1;
        auto c2 = (c1+1 < vecsize) ? c1+1 : c1;
        if(heaptree[c1] < heaptree[c2])
            max = c2;
        if(heaptree[max] > heaptree[k])
        {
            SwapIvec(heaptree,k,max);
            k = max;
        }
        else
            break;
    }
}

void ShuffleVec(vector<int> &ivec)
{
    uniform_int_distribution<unsigned> u(0,ivec.size()-1);
    //default_random_engine e_for_seed; //Use a random generator to generate a random seed
    //default_random_engine e_for_index(e_for_seed()); //Use the random seed to generate another random generator for randomizing
    //array index
    default_random_engine e_for_index;
    
    for(vector<int>::size_type i = 0; i < ivec.size(); ++i)
    {
        //Swap index i with random index u(e)
        int randindex1 = u(e_for_index);
        int randindex2 = u(e_for_index);
        //Sanity chk
        if((randindex1 < ivec.size()) && (randindex2 < ivec.size()))
            SwapIvec(ivec,randindex1,randindex2);
        //PrintSortVec(ivec);
    }
}

void MyQuickSort(vector<int> &ivec, vector<int>::size_type lo, vector<int>::size_type hi)
{
    static int reccnt = 0;
    cout << "Sort lo: " << lo << " hi: " << hi << " cnt: " << ++reccnt << endl;
    if((lo>=hi) || (lo>=ivec.size()) || (hi>=ivec.size()))
        return;
    auto j = MyQuickPartition(ivec,lo,hi);
    MyQuickSort(ivec,lo,j-1);
    MyQuickSort(ivec,j+1,hi);
    
    PrintSortVec(ivec);
}

vector<int>::size_type MyQuickPartition(vector<int> &ivec, vector<int>::size_type lo, vector<int>::size_type hi)
{
    static int reccnt = 0;
    cout << "                          Partition lo: " << lo << " hi: " << hi << " cnt: " << ++reccnt << endl;
    auto i = lo, j = hi+1;
    int v = ivec[lo]; //partition element
    
    while(true)
    {
        while((v > ivec[++i]) && (i <= hi)) {}
        while((v < ivec[--j]) && (j > lo)) {}
        if(i>=j)
            break;
        SwapIvec(ivec,i,j);
    }
    SwapIvec(ivec,lo,j);
    
    return j;
}

void SwapIvec(vector<int> &ivec, vector<int>::size_type lo, vector<int>::size_type hi)
{
    if(lo==hi)
        return;
    int temp = ivec[lo];
    ivec[lo] = ivec[hi];
    ivec[hi] = temp;
}

void MyMergeSort(vector<int> &ivec, vector<int>::size_type lo, vector<int>::size_type hi)
{
    static int reccnt = 0;
    cout << "Sort lo: " << lo << " hi: " << hi << " cnt: " << ++reccnt << endl;
    if(lo>=hi)
        return;
    auto mid = (lo+hi)/2;
    MyMergeSort(ivec,lo,mid);
    MyMergeSort(ivec,mid+1,hi);
    MyMerge(ivec,lo,hi);
    
    PrintSortVec(ivec);
}

void MyMerge(vector<int> &ivec, vector<int>::size_type lo, vector<int>::size_type hi)
{
    static int reccnt = 0;
    cout << "                          Merge lo: " << lo << " hi: " << hi << " cnt: " << ++reccnt << endl;
    
    if(ivec.size()<2)
        return;
    auto mid = (lo+hi)/2;
    
    vector<int> tempvec(ivec);
    
    decltype(lo) i = lo; //loop ctr for 1st half of the array. i.e lo to mid
    decltype(lo) j = mid+1; //loop ctr for 2nd half of the array. i.e mid+1 to hi
    
    for(decltype(lo) k = lo; k <= hi; ++k)
        if(j>hi)
            ivec[k] = tempvec[i++];
        else if(i>mid)
            ivec[k] = tempvec[j++];
        else if(tempvec[i] < tempvec[j])
            ivec[k] = tempvec[i++];
        else
            ivec[k] = tempvec[j++];
}

/*
 bool operator< (const string &s1, const string &s2)
 {
 if(s1.size()<s2.size())
 return true;
 else if(s1.size()>s2.size())
 return false;
 else
 {
 for(int i = 0; i < s1.size(); ++i)
 if (s1[i] < s2[i])
 return true;
 else if (s1[i] > s2[i])
 return false;
 return false;
 }
 }
 */

void PrintSetlist(const list<unordered_set<int> > &setlist)
{
    int setcnt = 0;
    for(const auto &i : setlist)
    {
        ++setcnt;
        cout << "Set#" << setcnt << " ";
        for(const auto &j : static_cast<unordered_set<int>>(i))
        {
            cout << j << " ";
        }
        cout << endl;
    }
}

void PrintSortVec(const vector<int> &ivec)
{
    for(const auto &i : ivec)
        cout << i << " ";
    cout << endl;
}

template <typename T>
bool Popstack(stack<T> &S, T &Val)
{
    if(!S.empty())
    {
        Val = S.top();
        S.pop();
        return true;
    }
    return false;
}

long long int Factorial(const int Temp)
{
    static int cnt = 0;
    ++cnt;
    
    cout << "N: " << Temp << " Cnt: " << cnt << endl;
    
    if(Temp<2)
        return 1;
    else
        return (Temp*Factorial(Temp-1));
}

int BinSearchRecursive(const int S, const vector<int> &A, const vector<int>::size_type lo, const vector<int>::size_type hi)
{
    static int cnt = 0;
    ++cnt;
    
    decltype(lo) mid = (lo+hi)/2;
    if(lo>hi)
    {
        cout << "Case 1: " << lo << " " << mid << " " << hi << " " << cnt << endl;
        return -1;
    }
    if(A[mid] == S)
    {
        cout << "Case 2: " << lo << " " << mid << " " << hi << " " << cnt << endl;
        return static_cast<int>(mid);
    }
    else if(S < A[mid])
    {
        cout << "Case 3: " << lo << " " << mid << " " << hi << " " << cnt << endl;
        return BinSearchRecursive(S,A,lo,mid-1);
    }
    else
    {
        cout << "Case 4: " << lo << " " << mid << " " << hi << " " << cnt << endl;
        return BinSearchRecursive(S,A,mid+1,hi);
    }
}

void WriteToQueue(queue<string> &WordQ, const string &Word)
{
    if( WordQ.size() >= 6 )
        WordQ.pop();
    WordQ.push(Word);
}

int f(void)
{
    cout << "f" << endl;
    return 1;
}
int g(void)
{
    cout << "g" << endl;
    return 1;
}
int h(void)
{
    cout << "h" << endl;
    return 1;
}
int j(void)
{
    cout << "j" << endl;
    return 1;
}

bool BinSearch(const int i, vector<int> *pvec)
{
    auto low = pvec->begin();
    auto hi = pvec->end();
    auto mid = pvec->begin() + pvec->size()/2;
    
    while( low != hi )
    {
        if(i == *mid)
            return true;
        else if (i < *mid)
            hi = mid;
        else
            low = mid + 1;
        mid = low + (hi - low)/2;
    }
    return false;
}
