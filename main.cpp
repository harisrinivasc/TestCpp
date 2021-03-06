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
#include <array>
//#include <cmath>
//#include "Sales_item.h"


//using std::string;
using namespace std;

//int reused = 42; // reused has global scope
typedef unsigned int uint;

class LinkListNodeType;
class LinkListNodeType
{
public:
    int val;
    shared_ptr<LinkListNodeType> next;
    LinkListNodeType() : val(0), next(nullptr) {}
    LinkListNodeType(const int i) : val(i), next(nullptr) {}
    ~LinkListNodeType()
    {
        cout << "Destruct: Value: " << val << endl;
    }
};

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

//depth first search
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

//breadth first search
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
void InsertKeyToBst(const int key, const string value);
shared_ptr<BstType> InsertKeyToBst(shared_ptr<BstType> Node, const int key, const string value);
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
bool PairSortComp(const pair<int,int> i, const pair<int,int> j);
bool PairSortComp2(const pair<int,int> i, const pair<int,int> j);
void IslandsBfs(const int row, const int col, vector<vector<bool>> &Visited, const vector<vector<int>> &IslandsInput);
int LargestNonAdjSum(const vector<int> &InputVec, const int idx, unordered_map<int,int> &Memo);
void AddElemMedian(const int i, priority_queue<int, vector<int>, greater<int>> &MinHeap, priority_queue<int> &MaxHeap);
void BalanceMedian(priority_queue<int, vector<int>, greater<int>> &MinHeap, priority_queue<int> &MaxHeap);
float GetMedian(const priority_queue<int, vector<int>, greater<int>> &MinHeap, const priority_queue<int> &MaxHeap);
void StringPerm(const vector<string> &Input, const string::size_type idx, const string &input_str, unordered_set<string> &StrPermOutputMap);
int GetNonuniformRandIdx(const int R, const vector<int> &CumProb, int low, int hi);
bool IsBst(const shared_ptr<BstType> Curr, const int Min, const int Max);
void InsertKeyToAnyTree(const int key, const string value);
shared_ptr<BstType> InsertKeyToAnyTree(shared_ptr<BstType> Node, const int key, const string value);
pair<bool,int> FindLargestBst(const shared_ptr<BstType> Curr, const int Min, const int Max);
pair<shared_ptr<BstType>,int> FindDeepestNodeTree(const shared_ptr<BstType> Curr, const int level);
void InvertBinaryTreeRecur(shared_ptr<BstType> Curr);
void InvertBinaryTreeIter(shared_ptr<BstType> Curr);
int HeightBinTree(const shared_ptr<BstType> Curr);
int SumOfAllNodes(const shared_ptr<BstType> Curr);
int MaxSumAnyLeafToRoot(const shared_ptr<BstType> Curr);
pair<shared_ptr<BstType>,int> SumLeafToRoot(const shared_ptr<BstType> Curr, const int key);
shared_ptr<BstType> LcaBt(const shared_ptr<BstType> Curr, const int A, const int B);
shared_ptr<BstType> LcaBst(const shared_ptr<BstType> Curr, const int A, const int B);
bool MyCompDecreasingOrder(const int i, const int j);
void AllPermOfNumVec(vector<vector<int>> &AllPermVecOp, vector<int> DigitCnt, int level, vector<int> TempVecOp);
bool HopFunction(const int curridx, const vector<int> &HopInput, unordered_map<int,bool> &HopMemo);
void ReverseWord(string &word, int begin, int end);
void MergeSortedVecs(const vector<int> &A, const vector<int> &B, vector<int> &Op);
pair<int,int> TreeLevelMinSum(const shared_ptr<BstType> Curr, int sum, int level);
bool IsSubtreeCheck(const shared_ptr<BstType> CurrT, const shared_ptr<BstType> CurrS, const shared_ptr<BstType> &BstHeadS );
void GenerateGrayCode(const int num, const int Len, unordered_set<int> &GrayCodeSet);
int MinSubsetFn(int Asum, int Bsum, int idx, const vector<int> &Input);
int MinPalinPartCnt(const string &Input, unordered_map<string,int> &Memo, const int i, const int j);
bool IsPalindrome(const string &Input, const int i, const int j);
int LenLongestDistint(const vector<int> &Input);
int NumSqInts(const int SqSumNum, const vector<int> &NumArray, int NumInts, int MaxIdx);
int MaxSumSubarray(const vector<int> &Input, const int start, const int end);
bool IsArrayEnd(const vector<int> &Input, const int curridx, unordered_map<int, bool> &Memo);
int ComputeInversions(const vector<pair<int,int>> &LineSegments);
int MergeSortInvCnt(vector<int> &Input, vector<int> &TempMerge, const int LeftStart, const int RightEnd);
int MergeSortInvCntCombine(vector<int> &Input, vector<int> &TempMerge, const int LeftStart, const int RightEnd);
int FindIdxSortedArray(const vector<vector<int>> &Input, const int FindNum, const int row, const int left, const int right);
int KthSmallest2d(const vector<vector<int>> &arr, const int k);
int SumOfNodesRepetitions(const shared_ptr<BstType> Curr, unordered_map<int, int> &SumCntMap);
void ReverseArray(vector<int> &Input, const int start, const int end);
bool TriangleErrorCheck(const vector<vector<int>> &Input);
pair<int,int> TriangleGetRight(const pair<int,int> Curr, const int M);
pair<int,int> TriangleGetLeft(const pair<int,int> Curr, const int M);
int TriangleMaxSum(const vector<vector<int>> &Input, const pair<int,int> Curr);
int FindMinElem(const vector<int> &Input, const int start, const int end);
int GetDigit(const int &Input, const int idx);
void PrintLinkedList(const shared_ptr<LinkListNodeType> Head);
shared_ptr<LinkListNodeType> AddNodeToLlHead(shared_ptr<LinkListNodeType> Head, shared_ptr<LinkListNodeType> Curr);
shared_ptr<LinkListNodeType> AddNodeToLlTail(shared_ptr<LinkListNodeType> Tail, shared_ptr<LinkListNodeType> Curr);
int FindCollatzCnt(const int N, unordered_map<int,int> &CollatzMap);
typedef pair<int,int> INT_PAIR;
struct CompareSortedKlists
{
    bool operator()(const INT_PAIR &a, const INT_PAIR &b)
    {
        return a.first < b.first;
    }
};
void FindIpAddr(const string input, const int stringidx, int IpByteVal, vector<int> &IpAddr, int IpBytePos);
void PrintDebugIpAddrData(const int &stringidx, const int &IpByteVal, const vector<int> &IpAddr, const int &IpBytePos, const int &LoopCtr);
shared_ptr<BstType> IsNodePresentBst(const shared_ptr<BstType> Curr, const int A);
int BstLevelCnt(const shared_ptr<BstType> Curr, const int A);
int FindContOnesLocFromLsb(int A);
void MultiplyMatrix(const vector<vector<int>> &A, const vector<vector<int>> &B, vector<vector<int>> &Output);
void ComputeMatrixPower(const vector<vector<int>> &A, const int N, vector<vector<int>> &Output);
int MultByPowerOf10(int num, int pow);
int FindKthElemToDelete(const int N, const int K);
string SplitSentance(const string &inString, const unordered_set<string> &WordDict, int idx);
void RunDfsArray(const int idx, vector<bool> &Visited, stack<char> &TopoSortOp, const vector<char> &CharArray);
bool CrackSafeDfs(const int &K, const int &N, unordered_set<string> &Visited, string &OpPwd);
int GetNumberOfDigits(int A);
bool MyCompSortForDcp228(const int A, const int B);
int CntNumAttempsEggs(const int numEggs, const int numFloors, unordered_map<string,int> &MemoEgg);
void FindAllKeypadComb(const int idx, const string &InNum, const vector<string> &NumAlphTable, vector<string> &KeypadOp, string &TempStr);
bool RearrangeRepeatChar(list<char> &CharList, string &OpStr);
void FindDecodeNumWays(const int Num, int &NumWays);
int FindQuotient(const int dividend, const int divisor);
int EarliestReachableStep(const vector<int> &MaxReach, const int low, const int high, const int find_idx);
bool MyStringSortComp(const char &A, const char &B);
bool FixedPntArray(const vector<int> &In, const int begin, const int end);
bool QuxRgb(list<int> &RgbIn);
void CheckCustDrinksComb(int currdrink, unordered_set<int> &CustSet, vector<unordered_set<int>> &MinDrinks, const vector<vector<int>> &DrinksList, unordered_set<int> &BestDrink);
bool MyCandSort(const pair<int,int> &A, const pair<int,int> &B);
int FindPosSortedArray(const vector<int> &In, const int begin, const int end, const int num);
vector<int> GetLockNextComb(const int curr, const int NumDigits);

class TempCl;
class TempCl
{
private:
    int x;
public:
    TempCl() : x(9)
    {
        cout << "Class constructed!" << endl;
    }
    ~TempCl()
    {
        cout << "Class destroyed!" << endl;
    }
    TempCl(const TempCl &other)
    {
        x = other.x;
        cout << "Class copied " << x << endl;
    }
};
TempCl CreateClass(void)
{
    TempCl *ptr = new TempCl;
    return *ptr;
}

class TrieType;
class TrieType : public enable_shared_from_this<TrieType>
{
private:
    bool IsEnd;
    unordered_map<char, shared_ptr<TrieType>> Child;
    void PrintAllWords(shared_ptr<TrieType> Curr, string &word_op)
    {
        if(Curr->IsEnd)
            cout << word_op << endl;
        for(auto itr = Curr->Child.cbegin(); itr != Curr->Child.cend(); ++itr)
        {
            word_op.push_back(itr->first);
            PrintAllWords(itr->second, word_op);
            word_op.pop_back();
        }
    }
public:
    TrieType()
    {
        IsEnd = false;
        Child.clear();
    }
    void AddWord(const string word_input)
    {
        shared_ptr<TrieType> Curr = shared_from_this(), TempNew;
        for(const auto &i:word_input)
        {
            if(Curr->Child.find(i) == Curr->Child.end())
            {
                TempNew = make_shared<TrieType>();
                Curr->Child.insert(make_pair(i, TempNew));
            }
            Curr = Curr->Child[i];
        }
        Curr->IsEnd = true;
    }
    void PrintAllWords(void)
    {
        string Output_word;
        PrintAllWords( shared_from_this(), Output_word );
    }
};

typedef struct
{
    int x;
    int y;
} RectCoordinateType;
istream & operator >> (istream &in, RectCoordinateType &r)
{
    in >> r.x >> r.y;
    return in;
}
ostream & operator << (ostream &out, const RectCoordinateType &r)
{
    out << "Xcoord " << r.x << " Ycoord " << r.y << endl;
    return out;
}

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
             InsertKeyToBst(key, value);
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
    //DCP168 - Given an N by N matrix, rotate it by 90 degrees clockwise.
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
    //Find maximum sum of contiguous subarray
    const vector<int> MaxSubarrayIn = {1,-3,4,-1,2,1};
    int MaxVal = 0, temp;
    for(auto i = 0; i < MaxSubarrayIn.size(); ++i)
        temp = FindArraySum(MaxSubarrayIn, i, MaxSubarrayIn.size()-1, MaxVal);
    cout << MaxVal << endl;
     */
    
    /*
    //DCP190 - Find maximum sum of contiguous subarray of a circular array
    //For example, given [8, -1, 3, 4], return 15 as we choose the numbers 3, 4, and 8 where the 8 is obtained from wrapping around.
    //Given [-4, 5, 1, 0], return 6 as we choose the numbers 5 and 1.
    vector<int> MaxSubarrayIn = {1,-3,4,-1,2,1};
    int MaxSumNoWrap = MaxSumSubarray(MaxSubarrayIn, 0, MaxSubarrayIn.size()-1);
    //To find sum of array with wraparound, compute sum and subtract MaxSum of negative subarray
    int TotalSum = 0;
    for(int i = 0; i < MaxSubarrayIn.size(); ++i)
    {
        TotalSum += MaxSubarrayIn[i];
        MaxSubarrayIn[i] = -MaxSubarrayIn[i];
    }
    int MaxSumNegArray = MaxSumSubarray(MaxSubarrayIn, 0, MaxSubarrayIn.size()-1);
    cout << max(MaxSumNoWrap, TotalSum-(-MaxSumNegArray)) << endl;
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
    //DCP158 - You are given an N by M matrix of 0s and 1s. Starting from the top left corner, how many ways are there to reach the bottom right corner?
    //You can only move right and down. 0 represents an empty space while 1 represents a wall you cannot walk through.
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
     //Given a list of possibly overlapping intervals, return a new list of intervals where all overlapping intervals have been merged.
     //The input list is not necessarily ordered in any way.
     //For example, given [(1, 3), (5, 8), (4, 10), (20, 25)], you should return [(1, 3), (4, 10), (20, 25)].
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
    
    /*
    //Generate random number from a given list and given distribution
    const vector<int> RandNums = {10,20,30,40,50};
    const vector<int> Prob = {3,2,1,1,3}; //probability distribution over 10
    vector<int> CumProb; //Cummulative probability
    
    CumProb.push_back(Prob[0]);
    for(int i = 1; i < Prob.size(); ++i)
        CumProb.push_back(CumProb[i-1]+Prob[i]);
    
    const int MaxProb = CumProb.back();
    vector<int> RandInxStats (static_cast<int>(Prob.size()), 0);
    
    for(int i = 0; i < 1000; ++i)
    {
        int R = (rand() % MaxProb)+1; //generate a random number in the range 1 to MaxProb
        ++RandInxStats[ GetNonuniformRandIdx(R,CumProb,0,CumProb.size()-1) ];
    }
    */
    
    /*
    //Check for number of braces to be removed to make it valid - DCP 86
    const string BraceIn = ")(((()))())(";
    stack<char> BraceStack;

    for(int i = 0; i < BraceIn.size(); ++i)
    {
        if(BraceIn[i] == '(')
            BraceStack.push(BraceIn[i]);
        else if(BraceIn[i] == ')')
        {
            if( !BraceStack.empty() && BraceStack.top() == '(' )
                BraceStack.pop();
            else
                BraceStack.push(BraceIn[i]);
        }
    }
    cout << BraceStack.size() << endl;
    */
    
    /*
    //Check if a tree is BST
    const vector<int> BstInputInt = {8,6,10,3,7};
    for(const auto &i : BstInputInt)
        InsertKeyToAnyTree(i, " ");
    PrintBst();
    cout << "Is BST " << IsBst(BstHead, INT_MIN, INT_MAX) << endl;
    */
    
    /*
    //Find the largest BST in a given tree/subtree - DCP #93
    const vector<int> TreeInputInt = {8,6,10,3,7};
    for(const auto &i : TreeInputInt)
        InsertKeyToAnyTree(i, " ");
    PrintBst();
    pair<bool,int> ret = FindLargestBst(BstHead, INT_MIN, INT_MAX);
    cout << "LargestBST " << ret.first << " " << ret.second << endl;
    */
    
    /*
    //Find the deepest node in tree - DCP 80
    const vector<int> TreeInputInt = {1,2,3,4,5,6,7};
    for(const auto &i : TreeInputInt)
        InsertKeyToBst(i, " ");
    PrintBst();
    pair<shared_ptr<BstType>,int> ret = FindDeepestNodeTree(BstHead, 1);
    cout << "Deepest node " << (ret.first ? ret.first->key : -1) << " level " << ret.second << endl;
    */
    
    /*
    //Invert binary tree - DCP 83
    const vector<int> TreeInputInt = {1,2,3,4,5,6,7};
    for(const auto &i : TreeInputInt)
        InsertKeyToBst(i, " ");
    PrintBst();
    InvertBinaryTreeRecur(BstHead);
    PrintBst();
    InvertBinaryTreeIter(BstHead);
    PrintBst();
    */
    
    /*
    //Find the max height of a tree
    const vector<int> TreeInputInt = {1,2,3,4,5,6,7,8};
    for(const auto &i : TreeInputInt)
        InsertKeyToAnyTree(i, " ");
    PrintBst();
    cout << "Tree height " << HeightBinTree(BstHead) << endl;
    */
    
    //Find least common ancestor for binary tree
    
    //Find least common ancestor for binary search tree
    
    /*
    //Find sum from a given leaf to root - tree
    const vector<int> TreeInputInt = {1,2,3,4,5,6,7};
    for(const auto &i : TreeInputInt)
        InsertKeyToAnyTree(i, " ");
    PrintBst();
    cout << "sum from leaf to root " << SumLeafToRoot(BstHead,4).second << endl;
    */
    
    /*
    //Find max. sum of path from any leaf to root - tree
    //Find sum of all elements under a root tree
    const vector<int> TreeInputInt = {1,2,3,4,5,6,7};
    for(const auto &i : TreeInputInt)
        InsertKeyToAnyTree(i, " ");
    PrintBst();
    cout << "Sum of all nodes " << SumOfAllNodes(BstHead) << endl;
    cout << "Max sum from any leaf to root " << MaxSumAnyLeafToRoot(BstHead) << endl;
    */
    
    /*
    //Lowest common ancestor (LCA) of Binary tree
    //Assumptions: both nodes are present in BT, no repetitions
    const vector<int> TreeInputInt = {1,2,3,4,5,6,7};
    for(const auto &i : TreeInputInt)
        InsertKeyToAnyTree(i, " ");
    PrintBst();
    shared_ptr<BstType> Ret = LcaBt(BstHead, 7,4);
    cout << "BT BST " << (Ret ? Ret->key : -1) << endl;
    */
    
    /*
    //Lowest common ancestor (LCA) of Binary search tree
    //Assumptions: both nodes are present in BT, no repetitions
    const vector<int> TreeInputInt = {4,2,1,3,6,5,7};
    for(const auto &i : TreeInputInt)
        InsertKeyToBst(i, " ");
    PrintBst();
    shared_ptr<BstType> Ret = LcaBst(BstHead, 7,13);
    cout << "BT BST " << (Ret ? Ret->key : -1) << endl;
    */
    
    /*
    //Find next greater permutation of a given number - DCP 95/DCP 205
    const int NextPermNum = 291874;
    //Convert given number to vector of digits
    vector<int> NumToVec;
    int temp = NextPermNum;
    while(temp)
    {
        NumToVec.push_back(temp % 10);
        temp /= 10;
    }
    //Find first digit that is < prev digit
    int low_idx = -1, high_idx = -1;
    for(int i = 1; i < NumToVec.size(); ++i)
        if(NumToVec[i] < NumToVec[i-1])
        {
            low_idx = i;
            break;
        }
    if(-1 == low_idx)
    {
        //current number is the max of all permutations.
        //return the smallest permutation
        low_idx = 0;
        high_idx = NumToVec.size()-1;
    }
    else
    {
        //Identify the smallest digit between low_idx-1 to 0
        high_idx = low_idx-1;
        for(int i = low_idx-2; i >= 0; --i)
            if(NumToVec[i] <= NumToVec[i+1])
                high_idx = i;
    }
    //swap low and high idx digits
    temp = NumToVec[low_idx];
    NumToVec[low_idx] = NumToVec[high_idx];
    NumToVec[high_idx] = temp;
    //Sort all digits between 0 to low_idx-1 in descending order
    sort(NumToVec.begin(), NumToVec.begin()+low_idx, MyCompDecreasingOrder);
        
    //Convert vector of digits to number
    int NumOutput = 0;
    for(int i = 0; i < NumToVec.size(); ++i)
        NumOutput += (NumToVec[i] * pow(10,i));
    cout << "Next perm " << NumOutput << endl;
    */
    
    /*
    //all permutations of a string - DCP 96
    const int AllPermNum = 11874;
    //Convert given number to vector of digits
    vector<int> NumToVec;
    int temp = AllPermNum;
    while(temp)
    {
        NumToVec.push_back(temp % 10);
        temp /= 10;
    }
    //Cnt the num. of occurence of each digit. Use a generic 128 ASCII vector
    vector<int> DigitCnt (128,0);
    for(const auto &i : NumToVec)
        ++DigitCnt['0'+i];
    vector<vector<int>> AllPermVecOp;
    vector<int> TempVecOp;
    AllPermOfNumVec(AllPermVecOp,DigitCnt,NumToVec.size()-1,TempVecOp);
    //Convert vector of digits to number
    for(const auto &i : AllPermVecOp)
    {
        int NumOutput = 0;
        for(int j = 0; j < i.size(); ++j)
            NumOutput += (i[i.size()-j-1] * pow(10,j));
        cout << NumOutput << endl;
    }
    */
    
    /*
    //Length of longest increasing consecutive elements sequence in an array - DCP 99
    //O(nlogn) complexity using heap
    const vector<int> LICSinput = {100,4,200,1,3,2,201,101,102,202,89,99,204,206,203,205};
    priority_queue<int> licsQ (LICSinput.begin(), LICSinput.end());
    priority_queue<int> MaxLenQ;
    int Len = 1;
    int prev = licsQ.top();
    licsQ.pop();
    while(!licsQ.empty())
    {
        int curr = licsQ.top();
        licsQ.pop();
        if(prev == curr+1)
            ++Len;
        else
        {
            MaxLenQ.push(Len);
            Len = 1;
        }
        prev = curr;
    }
    MaxLenQ.push(Len);
    cout << "Max len " << MaxLenQ.top() << endl;
    */
    
    /*
    //Length of longest increasing consecutive elements sequence in an array - DCP 99
    //O(n) complexity using set
    const vector<int> LICSinput = {100,4,200,1,3,2,201,101,102,202,89,99,204,206,203,205,5,199};
    unordered_set<int> licsSET (LICSinput.begin(), LICSinput.end());
    int MaxLen = INT_MIN;
    int array_idx = 0;
    while(!licsSET.empty() && array_idx < LICSinput.size())
    {
        int CurrLen = 1;
        int curr_elem = LICSinput[array_idx];
        int left = curr_elem-1;
        int right = curr_elem+1;
        ++array_idx;
        
        while(licsSET.find(left) != licsSET.end())
        {
            ++CurrLen;
            licsSET.erase(left);
            --left;
        }
        while(licsSET.find(right) != licsSET.end())
        {
            ++CurrLen;
            licsSET.erase(right);
            ++right;
        }
        
        if(MaxLen < CurrLen)
            MaxLen = CurrLen;
        licsSET.erase(curr_elem);
    }
    cout << "Max len " << MaxLen << endl;
     */
    
    /*
    //Least number of steps between given sequence of co-ordinates - DCP 100
    const vector<pair<int,int>> SeqSteps = {make_pair(0,0), make_pair(1,7), make_pair(2,0), make_pair(4,8), make_pair(0,1)};
    int stepcnt = 0;
    int xdiff, ydiff;
    for(int i = 1; i < SeqSteps.size(); ++i)
    {
        xdiff = abs(SeqSteps[i].first - SeqSteps[i-1].first);
        ydiff = abs(SeqSteps[i].second - SeqSteps[i-1].second);
        stepcnt += min(xdiff,ydiff);
        stepcnt += abs(xdiff-ydiff);
    }
    cout << stepcnt << endl;
    */
    
    /*
    //Find 2 prime integers that sum up to a given even number - DCP 101
    vector<int> PrimeNums = {2}; //look-up table for all prime numbers between 2 and PrimeNumMax
    const int PrimeNumMax = 100000;
    int IsPrime = true;
    for(unsigned int i = 3; i < PrimeNumMax; i=i+2)
    {
        IsPrime = true;
        int SqRoot = sqrt(i);
        for(const auto &j : PrimeNums)
        {
            if(0 == i%j && j <= SqRoot)
            {
                IsPrime = false;
                break;
            }
        }
        if(IsPrime)
            PrimeNums.push_back(i);
    }
    const unordered_set<int> PrimeSet (PrimeNums.begin(), PrimeNums.end());
    
    const int PrimeSumInput = 787976;
    int diff = 0;
    if(PrimeSumInput > PrimeNumMax)
        cerr << "Prime sum input > max range" << endl;
    else
    {
        for(const auto &i : PrimeNums)
        {
            diff = PrimeSumInput-i;
            if(PrimeSet.find(diff) != PrimeSet.end())
            {
                cout << i << " " << diff << endl;
                break;
            }
        }
    }
    */
    
    /*
    //Find contiguous subarray that sum to a given number - DPC  102
    const vector<int> ContigSumInput = {22,0,1,144,6,20};
    const int FindSum = 44;
    int CurrSum = FindSum;
    int startidx = 0, endidx = 0;
    for(int i = 0; i < ContigSumInput.size(); ++i)
    {
        CurrSum -= ContigSumInput[i];
        endidx = i;

        if(CurrSum < 0)
        {
            while((CurrSum < 0) && (startidx <= endidx))
            {
                CurrSum += ContigSumInput[startidx];
                startidx++;
            }
        }
        if(0 == CurrSum)
        {
            cout << "Subarray found!! " << startidx << " " << endidx << endl;
            break;
        }
    }
    if(CurrSum != 0)
        cout << "Contiguous subarray not found" << endl;
    */
    
    /*
    //hop steps - DCP 106
    const vector<int> HopInput = {3,1,0,1};
    unordered_map<int,bool> HopMemo;
    cout << HopFunction(0, HopInput, HopMemo) << endl;
    */
    
    /*
    //reverse a line - DCP 113
    string ReverseLine = "hello  world here";
    ReverseWord(ReverseLine, 0, ReverseLine.size()-1); //reverse entire line
    cout << ReverseLine << endl;
    
    int begin = 0, end = 0;
    while(end < ReverseLine.size())
    {
        //Reverse word in-between spaces
        ++end;
        if((ReverseLine[end] == ' ') || (end == ReverseLine.size()))
        {
            ReverseWord(ReverseLine, begin, end-1);
            begin = end+1;
        }
    }
    cout << ReverseLine << endl;
    */
    
    /*
    //reverse a line but keep the order of delimiters - DCP 114
    const string ReverseLine = "hello/*world?here";
    string ReverseLineOp (ReverseLine.size(), ' ');
    queue<string> ReverseQ;
    stack<string> ReverseStack;
    const unordered_set<char> Delim = {'/', '*', '?', ' '};
    string TempS;
    
    int begin = 0, end = 0;
    while(end < ReverseLine.size())
    {
        if(Delim.find(ReverseLine[end]) != Delim.end())
        {
            TempS = ReverseLine.substr(ReverseLine.begin()+begin, ReverseLine.begin()+end-1);
            ReverseStack.push(TempS);
     }
     }
    */
    
    /*
    //find the smallest set of numbers that cover all the intervals - DPC 119
    vector<pair<int,int>> Intervals {make_pair(1, 3), make_pair(5, 8), make_pair(4, 10), make_pair(20, 25), make_pair(1,10)};
    vector<int> OutputIntervals = {Intervals[0].second};
    
    //Sort the input based on start times (first element in pair)
    sort(Intervals.begin(), Intervals.end(), PairSortComp2);
    
    for(int i = 1; i < Intervals.size(); ++i)
    {
        if(Intervals[i].first <= OutputIntervals.back())
            continue;
        if(OutputIntervals.size() > 1)
            OutputIntervals.pop_back();
        OutputIntervals.push_back(Intervals[i].first);
    }
     */
    
    /*
    //Square the sorted elements and give output in sorted order - DCP 118
    const vector<int> SquareIn = {-10,-8,-8,-7,-3,-1,0,7,8,8,10,11,12,13,17};
    int Negvecsize = 0;
    for(const auto &i: SquareIn)
        if(i < 0)
            ++Negvecsize;
        else
            break;
    vector<int> NegVec (Negvecsize,0), PosVec (SquareIn.size()-Negvecsize,0);
    vector<int> SquareOp;
    //Split input array into +ve and -ve vectors
    for(int i = 0; i < Negvecsize; ++i)
        NegVec[i] = abs(SquareIn[Negvecsize-i-1]);
    for(int i = Negvecsize; i < SquareIn.size(); ++i)
        PosVec[i-Negvecsize] = SquareIn[i];
    //Merge the +ve and -ve sorted vectors
    MergeSortedVecs(NegVec, PosVec, SquareOp);
    //Square the elements
    for(int i = 0; i < SquareOp.size(); ++i)
        SquareOp[i] *= SquareOp[i];
    */
    
    /*
    //Find level of min. sum of a tree - DCP 117
    const vector<int> TreeInputInt = {5,4,3,2,-10,6};
    for(const auto &i : TreeInputInt)
        InsertKeyToBst(i, " ");
    PrintBst();
    auto ret = TreeLevelMinSum(BstHead, 0, 0);
    cout << "Sum: " << ret.first << " Level: " << ret.second << endl;
    */
    
    /*
    //check if a tree 't' or its subtree has a same value and strcuture as another tree 's' - DCP 115
    shared_ptr<BstType> BstHeadT = nullptr, BstHeadS = nullptr;
    const vector<int> TreeInputT = {5,4,3,2,-10,7,6,8};
    for(const auto &i : TreeInputT)
        BstHeadT = InsertKeyToBst(BstHeadT,i, " ");
    PrintBst(BstHeadT);
    const vector<int> TreeInputS = {7,6,8};
    for(const auto &i : TreeInputS)
        BstHeadS = InsertKeyToBst(BstHeadS,i, " ");
    PrintBst(BstHeadS);
    cout << IsSubtreeCheck(BstHeadT, BstHeadS, BstHeadS) << endl;
    */
    
    /*
    //DCP148 - generate gray code
    const int GrayCodeLen = 4;
    unordered_set<int> GrayCodeSet;
    GrayCodeSet.insert(0);
    GenerateGrayCode(0, GrayCodeLen, GrayCodeSet);
    */
    
    /*
    //DCP164 - You are given an array of length n + 1 whose elements belong to the set {1, 2, ..., n}. By the pigeonhole principle, there must be a duplicate. Find it in linear time and space.
    const vector<int> PigeonHoleIn = {1,1,4,3,2,5};
    int array_sum = 0, expected_sum = 0, n = PigeonHoleIn.size()-1;
    for(const auto &i: PigeonHoleIn)
        array_sum += i;
    expected_sum = n*(n+1)/2;
    cout << array_sum-expected_sum << endl;
    */
    
    /*
    //DCP185 - Given two rectangles on a 2D graph, return the area of their intersection. If the rectangles don't intersect, return 0.
    //Given left bottom and top right indicies of each rectangle.
    pair<int,int> Rec1LB = make_pair(1,1), Rec1TR = make_pair(4,3), Rec2LB = make_pair(1,1), Rec2TR = make_pair(5,3);
    if(( Rec2LB.first >= Rec1TR.first ) ||
       ( Rec2LB.second >= Rec1TR.second ))
        cout << "0" << endl;
    else
    {
        int x = abs( Rec2LB.first - Rec1TR.first);
        int y = abs( Rec2LB.second - Rec1TR.second);
        cout << "area " << x*y << endl;
    }
     */
    
    /*
    //DCP184 - Given n numbers, find the greatest common denominator between them.
    //For example, given the numbers [42, 56, 14], return 14.
    vector<int> GcdIn = {42,57,14};
    int p = GcdIn[0], q, r;
    for(int i = 1; i < GcdIn.size(); ++i)
    {
        q = GcdIn[i];
        while(0 != q)
        {
            r = p % q;
            p = q;
            q = r;
        }
    }
    cout << "gcd: " << p << endl;
    */
    
    /*
    //DCP186 - Given an array of positive integers, divide the array into two subsets such that the difference between the sum of the subsets is as small as possible.
    //For example, given [5, 10, 15, 20, 25], return the sets {10, 25} and {5, 15, 20}, which has a difference of 5, which is the smallest possible difference.
    vector<int> SubsetDiffIn = {5,4,3,2,1};
    int TotalSum = 0;
    for(const auto &i:SubsetDiffIn)
        TotalSum += i;

    cout << MinSubsetFn(TotalSum,0,SubsetDiffIn.size()-1,SubsetDiffIn) << endl;
    */
    
    /*
    //>> and << operator overloading
    vector<RectCoordinateType> RectCoordIn;
    RectCoordinateType Temp;
    for(int i = 0; i < 2; ++i)
    {
        cout << "Enter x and y coordinate #" << i << endl;
        cin >> Temp;
        RectCoordIn.push_back(Temp);
    }
    for(const auto &i : RectCoordIn)
        cout << i;
    */
    
    /*
    //DCP181 - Given a string, split it into as few strings as possible such that each string is a palindrome.
    //For example, given the input string racecarannakayak, return ["racecar", "anna", "kayak"].
    //Given the input string abc, return ["a", "b", "c"].
    //Similar to matrix chain multiplication
    const string PalinPartitionIn = "abcdcba";
    unordered_map<string,int> PalinPartitionMemo;
    cout << MinPalinPartCnt(PalinPartitionIn, PalinPartitionMemo, 0, PalinPartitionIn.size()-1) << endl;
    */
    
    /*
    //DCP189 - Given an array of elements, return the length of the longest subarray where all its elements are distinct.
    //For example, given the array [5, 1, 3, 5, 2, 3, 4, 1], return 5 as the longest subarray of distinct elements is [5, 2, 3, 4, 1].
    const vector<int> DistinctSubarray = {5, 1, 3, 5, 2, 3, 4, 4};
    cout << LenLongestDistint(DistinctSubarray) << endl;
    */
    
    /*
    //DCP154 - Implement a stack API using only a heap. A stack implements the following methods:
    //push(item), which adds an element to the stack
    //pop(), which removes and returns the most recently added element (or throws an error if there is nothing on the stack)
    //Recall that a heap has the following operations:
    //push(item), which adds a new key to the heap
    //pop(), which removes and returns the max value of the heap
    typedef pair<int,int> MyPair;
    class MyStackHeapComp
    {
    public:
        bool operator()(const MyPair &A, const MyPair &B)
        {
            return (A.first < B.first);
        }
    };
    
    class MyStack
    {
    private:
        priority_queue<MyPair, vector<MyPair>, MyStackHeapComp> MyHeap;
        //priority_queue<MyPair, vector<MyPair>, less<MyPair>> MyHeap;
        //priority_queue<MyPair> MyHeap;
    public:
        void push(const int a)
        {
            int CurrSize = MyHeap.size();
            MyHeap.push(make_pair(CurrSize,a));
        }
        void pop()
        {
            MyHeap.pop();
        }
        int top()
        {
            return MyHeap.top().second;
        }
    };
    MyStack temp_stack;
    temp_stack.push(2);
    temp_stack.push(3);
    temp_stack.push(0);
    temp_stack.pop();
    cout << temp_stack.top() << endl;
    */
    
    /*
    //DCP156 - Given a positive integer n, find the smallest number of squared integers which sum to n.
    //For example, given n = 13, return 2 since 13 = 3^2 + 2^2 = 9 + 4.
    //Given n = 27, return 3 since 27 = 3^2 + 3^2 + 3^2 = 9 + 9 + 9.
    const int SqSumNum = 3;
    vector<int> NumArray;
    for(int i = 1; i <= sqrt(SqSumNum); ++i)
        NumArray.push_back(i);
    
    cout << NumSqInts(SqSumNum, NumArray, 0, NumArray.size()-1) << endl;
     */
    
    /*
    //DCP161 - Given a 32-bit integer, return the number with its bits reversed.
    //For example, given the binary number 1111 0000 1111 0000 1111 0000 1111 0000, return 0000 1111 0000 1111 0000 1111 0000 1111.
    const unsigned int RevIn = INT_MAX;
    unsigned int RevOp = 0, InMask, OpMask;
    for(int i = 0; i < sizeof(unsigned int)*8; ++i)
    {
        InMask = 1 << i;
        if(RevIn & InMask)
        {
            OpMask = 1 << (sizeof(unsigned int)*8 - i - 1);
            RevOp |= OpMask;
        }
    }
    cout << RevIn << " " << RevOp << endl;
     */
    
    /*
    //DCP165 - Given an array of integers, return a new array where each element in the new array is the number of smaller elements to the right of that element in the original input array.
    //For example, given the array [3, 4, 9, 6, 1], return [1, 1, 2, 1, 0]
    const vector<int> SmallElemRightIn = {3, 4, 9, 6, 1};
    vector<int> NumElemOp (SmallElemRightIn.size(), 0);
    set<int> SmallElemSet;
    for(int i = SmallElemRightIn.size()-1; i >= 0; --i)
    {
        //Add element from the right to set (balanced binary search tree usin AVL)
        SmallElemSet.insert(SmallElemRightIn[i]);
        
        //find the index of the element in the BST
        auto itr = SmallElemSet.lower_bound(SmallElemRightIn[i]);
        
        //find the distance of this element from the start of the BST and store in output array
        NumElemOp[i] = distance(SmallElemSet.begin(), itr);
    }
    */
    
    /*
    //DCP180 - Given a stack of N elements, interleave the first half of the stack with the second half reversed using only one other queue. This should be done in-place.
    //Recall that you can only push or pop from a stack, and enqueue or dequeue from a queue.
    //For example, if the stack is [1, 2, 3, 4, 5], it should become [1, 5, 2, 4, 3]. If the stack is [1, 2, 3, 4], it should become [1, 4, 2, 3].
    const vector<int> InNum = {1, 2, 3, 4};
    stack<int> InterStack;
    for(const auto &i : InNum)
        InterStack.push(i);
    queue<int> InterQ;
    int startidx = 1, i;
    
    while(startidx < InNum.size())
    {
        i = startidx;
        while(i < InNum.size())
        {
            InterQ.push(InterStack.top());
            InterStack.pop();
            ++i;
        }
        while(!InterQ.empty())
        {
            InterStack.push(InterQ.front());
            InterQ.pop();
        }
        ++startidx;
    }
    */
    
    /*
    //DCP192 - You are given an array of nonnegative integers. Let's say you start at the beginning of the array and are trying to advance to the end. You can advance at most, the number of steps that you're currently on. Determine whether you can get to the end of the array.
    //For example, given the array [1, 3, 1, 2, 0, 1], we can go from indices 0 -> 1 -> 3 -> 5, so return true.
    //Given the array [1, 2, 1, 0, 0], we can't reach the end, so return false.
    const vector<int> EndArrayIn = {1, 3, 1, 1, 0, 1};
    unordered_map<int, bool> Memo;
    cout << IsArrayEnd(EndArrayIn, 0, Memo) << endl;
    */
    
    /*
    //Maximum difference between elements in an array
    const vector<int> MaxDiffInput = {2, 3, 10, 6, -4, 8, 1};
    int MinElem = MaxDiffInput[0], MaxDiff = MaxDiffInput[1] - MaxDiffInput[0];
    for(int i = 1; i < MaxDiffInput.size(); ++i)
    {
        MaxDiff = max(MaxDiff, MaxDiffInput[i]-MinElem);
        MinElem = min(MinElem, MaxDiffInput[i]);
    }
    cout << MaxDiff << " " << MinElem << endl;
    */
    
    /*
    //DCP193 - Given a array of numbers representing the stock prices of a company in chronological order, write a function that calculates the maximum profit you could have made from buying and selling that stock. You're also given a number fee that represents a transaction fee for each buy and sell transaction.
    //You must buy before you can sell the stock, but you can make as many transactions as you like.
    //For example, given [1, 3, 2, 8, 4, 10] and fee = 2, you should return 9, since you could buy the stock at 1 dollar, and sell at 8 dollars, and then buy it at 4 dollars and sell it at 10 dollars. Since we did two transactions, there is a 4 dollar fee, so we have 7 + 6 = 13 profit minus 4 dollars of fees.
    const vector<int> StockPrices = {1, 3, 2, 8, 4, 10};
    const int Fee = 2;
    int CashInHand = -StockPrices[0];
    int Profit = 0;
    for(int i = 1; i < StockPrices.size(); ++i)
    {
        Profit = max(Profit, StockPrices[i]+CashInHand-Fee);
        CashInHand = max(CashInHand, Profit-StockPrices[i]);
    }
    cout << Profit << endl;
    */
    
    /*
    //Given an array of stock prices, determine the maximum profit by buying and selling stock. Similar to DCP193 but without
    //transaction fee
    const vector<int> StockPrices = {10, 22, 5, 75, 65, 80};//{1, 3, 2, 8, 4, 10};
    vector<pair<int,int>> BuySell;
    int idx = 0, buy, sell;
    while (idx < StockPrices.size())
    {
        //Find local minima
        while( (StockPrices[idx+1] <= StockPrices[idx]) && (idx < StockPrices.size()-1) )
            ++idx;
        if(idx == StockPrices.size()-1)
            break;
        buy = idx++;
        
        //Find local maxima
        while( (StockPrices[idx] >= StockPrices[idx-1]) && (idx < StockPrices.size()) )
            ++idx;
        sell = idx-1;
        
        BuySell.push_back(make_pair(buy, sell));
    }
    */
    
    //Maximize profit by buying/selling atmost 'k' times
    //const int Ktimes = 2;

    /*
    //DCP194 - Suppose you are given two lists of n points, one list p1, p2, ..., pn on the line y = 0 and the other list q1, q2, ..., qn on the line y = 1. Imagine a set of n line segments connecting each point pi to qi. Write an algorithm to determine how many pairs of the line segments intersect.
    //pair.first indicates 'q' points on the line y = 1, pair.second indicates 'p' points on line y = 0
    vector<pair<int,int>> LineSegments {{make_pair(4,5)}, {make_pair(7,8)}, {make_pair(3,3)}, {make_pair(6,2)}, {make_pair(2,1)},
        {make_pair(1,7)}, {make_pair(4,7)}, {make_pair(10,12)}};
    
    //Sort based on y = 1 line
    sort(LineSegments.begin(), LineSegments.end(), PairSortComp);
    
    //Compute total number of inversions in the 'pair.second' elements using MergeSort. Total number of inversions indicate the total
    //number of intersections in the ling segment
    cout << ComputeInversions(LineSegments) << endl;
    */
    
    /*
    //DCP195 - Let A be an N by M matrix in which every row and every column is sorted.
    //Given i1, j1, i2, and j2, compute the number of elements of M smaller than M[i1, j1] and larger than M[i2, j2].
    const vector<vector<int>> Matrix2dIn = {{1, 3, 7, 10, 15, 20},
                                            {2, 6, 9, 14, 22, 25},
                                            {3, 8, 10, 15, 25, 30},
                                            {10, 11, 12, 23, 30, 35},
                                            {20, 25, 30, 35, 40, 45}};
    const pair<int,int> LessThan = make_pair(1,1), GreaterThan = make_pair(3,3);
    int TotalCnt = 0;
    
    //Cnt total elements less than given number
    for(int row = 0; row < Matrix2dIn.size(); ++row)
    {
        if(row >= LessThan.first)
            TotalCnt += (FindIdxSortedArray(Matrix2dIn, Matrix2dIn[LessThan.first][LessThan.second], row, 0, LessThan.second));
        else
            TotalCnt += (FindIdxSortedArray(Matrix2dIn, Matrix2dIn[LessThan.first][LessThan.second], row, LessThan.second, Matrix2dIn[row].size()-1));
    }
    
    //Cnt total elements greater than given number
    for(int row = 0; row < Matrix2dIn.size(); ++row)
    {
        if(row <= GreaterThan.first)
            TotalCnt += (Matrix2dIn[row].size() - FindIdxSortedArray(Matrix2dIn, Matrix2dIn[GreaterThan.first][GreaterThan.second], row, GreaterThan.second, Matrix2dIn[row].size()-1));
        else
            TotalCnt += (Matrix2dIn[row].size() - FindIdxSortedArray(Matrix2dIn, Matrix2dIn[GreaterThan.first][GreaterThan.second], row, 0, GreaterThan.second));
    }

    cout << TotalCnt << endl;
    */
    
    /*
    //compute kth smallest element in a sorted 2D array
    const vector<vector<int>> Matrix2dIn = {{1, 3, 4, 5, 15, 20},
                                            {2, 6, 9, 14, 22, 25},
                                            {3, 8, 10, 15, 25, 30},
                                            {10, 11, 12, 23, 30, 35},
                                            {20, 25, 30, 35, 40, 45}};
    cout << KthSmallest2d(Matrix2dIn, 2) << endl;
    */
    
    /*
    //DCP196 - Given the root of a binary tree, find the most frequent subtree sum. The subtree sum of a node is the sum of all values under a node, including the node itself.
    const vector<int> SumTree = {5,5,5,5,5,5,5};
    for(const auto &i : SumTree)
        InsertKeyToAnyTree(i, " ");
    PrintBst();
    unordered_map<int, int> SumCntMap;
    (void)SumOfNodesRepetitions(BstHead, SumCntMap);
    
    int MaxSum = 0, MaxSumCnt = 0;
    for(auto i = SumCntMap.cbegin(); i != SumCntMap.cend(); ++i)
    {
        if(MaxSumCnt < i->second)
        {
            MaxSumCnt = i->second;
            MaxSum = i->first;
        }
    }
    cout << MaxSum << " " << MaxSumCnt << endl;
     */
    
    /*
    //DCP197 - Given an array and a number k that's smaller than the length of the array, rotate the array to the right k elements in-place.
    vector<int> RotArray;
    for(int i = 0; i < 10; ++i)
        RotArray.push_back(i);
    const int RotK = 9;
    
    //Reverse entire array
    ReverseArray(RotArray, 0, RotArray.size()-1);
    
    //Reverse between 0 to K-1
    ReverseArray(RotArray, 0, RotK-1);
    
    //Reverse between k to size()-1
    ReverseArray(RotArray, RotK, RotArray.size()-1);
    
    for(const auto &i : RotArray)
        cout << i << " ";
    cout << endl;
    */
    
    /*
    //Rotate array by 1
    vector<int> RotArray;
    for(int i = 0; i < 10; ++i)
        RotArray.push_back(i);
    
    int temp = RotArray.back();
    for(int i = RotArray.size()-1; i > 0; --i)
        RotArray[i] = RotArray[i-1];
    RotArray[0] = temp;
    */
    
    /*
    //DCP198 - Given a set of distinct positive integers, find the largest subset such that every pair of elements in the subset (i, j) satisfies either i % j = 0 or j % i = 0.
    //For example, given the set [3, 5, 10, 20, 21], you should return [5, 10, 20]. Given [1, 3, 6, 24], return [1, 3, 6, 24].
    vector<int> Mod0Input = {1,2,3,5,7,11,13,17,19};//{12, 3, 5, 10, 20, 21, 6, 24};
    sort(Mod0Input.begin(), Mod0Input.end());
    vector<int> Mod0Output, TempVec;
    vector<bool> Visited(Mod0Input.size(), false);
    
    for(int i = Mod0Input.size()-1; i >= 0; --i)
    {
        TempVec.clear();
        TempVec.push_back(Mod0Input[i]);
        
        for(int j = i-1; j >= 0 && !Visited[j]; --j)
        {
            if(TempVec.back() % Mod0Input[j] == 0)
            {
                Visited[j] = true;
                TempVec.push_back(Mod0Input[j]);
            }
        }
        
        Visited[i] = true;
        if(TempVec.size() > Mod0Output.size())
            Mod0Output = TempVec;
    }
    */
    
    /*
    //DCP201 - You are given an array of arrays of integers, where each array corresponds to a row in a triangle of numbers. For example, [[1], [2, 3], [1, 5, 1]] represents the triangle:
    //We define a path in the triangle to start at the top and go down one row at a time to an adjacent value, eventually ending with an entry on the bottom row. For example, 1 -> 3 -> 5. The weight of the path is the sum of the entries.
    vector<vector<int>> TriangleInput;
    vector<int> Temp = {1};
    TriangleInput.push_back(Temp);
    Temp = {2,3};
    TriangleInput.push_back(Temp);
    Temp = {1,5,1};
    TriangleInput.push_back(Temp);
    
    if(TriangleErrorCheck(TriangleInput))
        cout << TriangleMaxSum(TriangleInput, make_pair(0,0)) << endl;
    else
        cout << "Error triangle input" << endl;
    */
    
    /*
    //DCP203 - Suppose an array sorted in ascending order is rotated at some pivot unknown to you beforehand. Find the minimum element in O(log N) time. You may assume the array does not contain duplicates.
    //For example, given [5, 7, 10, 3, 4], return 3.
    const vector<int> RotArrayInput = {5, 7, 10, 13, 14, 1, 2};
    cout << FindMinElem(RotArrayInput, 0, RotArrayInput.size()-1) << endl;
    */
    
    /*
    //DCP202 - Write a program that checks whether an integer is a palindrome. For example, 121 is a palindrome, as well as 888. 678 is not a palindrome. Do not convert the integer into a string.
    const int PalindromeNum = 543312345;
    int NumLen = 1;
    int quotient = PalindromeNum/10;
    bool IsNumPalin = true;
    
    while(quotient)
    {
        quotient /= 10;
        ++NumLen;
    }
    
    for(int i = 0; i <= NumLen/2; ++i)
        if( GetDigit(PalindromeNum,i) != GetDigit(PalindromeNum,NumLen-i-1) )
        {
            IsNumPalin = false;
            break;
        }
    cout << IsNumPalin << endl;
    */
    
    /*
    //DCP200 - Let X be a set of n intervals on the real line. We say that a set of points P "stabs" X if every interval in X contains at least one point in P. Compute the smallest set of points that stabs X.
    //For example, given the intervals [(1, 4), (4, 5), (7, 9), (9, 12)], you should return [4, 9].
    vector<pair<int,int>> IntervalsIn {{make_pair(8,9)}, {make_pair(2,6)}, {make_pair(1,8)}, {make_pair(1,5)}, {make_pair(5,7)}};
    //vector<pair<int,int>> IntervalsIn {{make_pair(9,12)}, {make_pair(4,5)}, {make_pair(7,9)}, {make_pair(1,4)}};
    
    //Sort based on .first
    sort(IntervalsIn.begin(), IntervalsIn.end(), PairSortComp2);
    
    vector<pair<int,int>> TempInterval;
    TempInterval.push_back(IntervalsIn[0]);
    for(int i = 1; i < IntervalsIn.size(); ++i)
        if(IntervalsIn[i].first != IntervalsIn[i-1].first)
            TempInterval.push_back(IntervalsIn[i]);
    
    vector<int> Ppoints;
    pair<int,int> MinMaxPoints = make_pair(INT_MAX,INT_MAX);
    
    for(int i = 0; i < TempInterval.size(); i++)
    {
        if(MinMaxPoints.first == INT_MAX)
            MinMaxPoints = TempInterval[i];
        else if(TempInterval[i].first == MinMaxPoints.second)
        {
            Ppoints.push_back(MinMaxPoints.second);
            MinMaxPoints = make_pair(INT_MAX,INT_MAX);
        }
        else if(TempInterval[i].first > MinMaxPoints.second)
        {
            Ppoints.push_back(MinMaxPoints.second);
            MinMaxPoints = TempInterval[i];
        }
        else
            MinMaxPoints = make_pair(max(TempInterval[i].first, TempInterval[i-1].first), min(TempInterval[i].second, TempInterval[i-1].second));
    }
    if(MinMaxPoints.first != INT_MAX)
        Ppoints.push_back(MinMaxPoints.first);
    */
    
    /*
    //DCP199 - Given a string of parentheses, find the balanced string that can be produced from it using the minimum number of insertions and deletions. If there are multiple solutions, return any of them.
    //For example, given "(()", you could return "(())". Given "))()(", you could return "()()()()"
    //Similar to DCP86
    const string BraceIn = "(()";//)(((()))())(";
    stack<pair<char,int>> TempBraceStack;
    
    for(int i = 0; i < BraceIn.size(); ++i)
    {
        if(BraceIn[i] == '(')
            TempBraceStack.push(make_pair(BraceIn[i],i));
        else if(BraceIn[i] == ')')
        {
            if( !TempBraceStack.empty() && TempBraceStack.top().first == '(' )
                TempBraceStack.pop();
            else
                TempBraceStack.push(make_pair(BraceIn[i],i));
        }
    }
    
    stack<char> OutputBraceStack;
    int idx = BraceIn.size()-1;
    while(idx >= 0)
    {
        if(!TempBraceStack.empty())
            if(TempBraceStack.top().second == idx)
            {
                //push closing brace first before opening brace
                OutputBraceStack.push(')');
                OutputBraceStack.push('(');
                TempBraceStack.pop();
                --idx;
                continue;
            }
        OutputBraceStack.push(BraceIn[idx--]);
    }
    
    string BraceOutput (OutputBraceStack.size(), ' ');
    idx = 0;
    while(idx < BraceOutput.size() && !OutputBraceStack.empty())
    {
        BraceOutput[idx++] = OutputBraceStack.top();
        OutputBraceStack.pop();
    }
    cout << BraceOutput << endl;
    */
    
    /*
    //DCP208 - Given a linked list of numbers and a pivot k, partition the linked list so that all nodes less than k come before nodes greater than or equal to k.
    //For example, given the linked list 5 -> 1 -> 8 -> 0 -> 3 and k = 3, the solution could be 1 -> 0 -> 5 -> 8 -> 3.
    const vector<int> LLinput {5,1,8,3,0,3,-1,6,9,4,-2};
    const int LLpivot = 3;
    shared_ptr<LinkListNodeType> LLHead = nullptr, TempLLnode = nullptr;
    
    //create Linked list by adding new elements to the head
    for(const auto &i:LLinput)
    {
        TempLLnode = make_shared<LinkListNodeType>(i);
        LLHead = AddNodeToLlHead(LLHead, TempLLnode);
    }
    PrintLinkedList(LLHead);
    
    //Split the list into lesser than list, and greater than or equal to list
    shared_ptr<LinkListNodeType> HeadGE = nullptr, TailGE = nullptr, HeadLess = nullptr, TailLess = nullptr;
    while(LLHead)
    {
        TempLLnode = LLHead;
        LLHead = LLHead->next;
        if(TempLLnode->val < LLpivot)
        {
            TailLess = AddNodeToLlTail(TailLess, TempLLnode);
            if(nullptr == HeadLess)
                HeadLess = TailLess;
        }
        else
        {
            TailGE = AddNodeToLlTail(TailGE, TempLLnode);
            if(nullptr == HeadGE)
                HeadGE = TailGE;
        }
    }
    PrintLinkedList(HeadLess);
    PrintLinkedList(HeadGE);
    
    //Combine the 2 lists
    LLHead = HeadLess;
    TailLess->next = HeadGE;
    PrintLinkedList(LLHead);
    */
    
    /*
    //DCP209 - Write a program that computes the length of the longest common subsequence of three given strings. For example, given "epidemiologist", "refrigeration", and "supercalifragilisticexpialodocious", it should return 5, since the longest common subsequence is "eieio".
    const string LcsIn1 = "epidemiologist", LcsIn2 = "refrigeration", LcsIn3 = "supercalifragilisticexpialodocious";
    vector<vector<vector<int>>> LcsDynArray (LcsIn1.size()+1, vector<vector<int>>(LcsIn2.size()+1, vector<int>(LcsIn3.size()+1,0)));
    for(int i = 1; i <= LcsIn1.size(); ++i)
        for(int j = 1; j <= LcsIn2.size(); ++j)
            for(int k = 1; k <= LcsIn3.size(); ++k)
                if(LcsIn1[i-1] == LcsIn2[j-1] && LcsIn2[j-1] == LcsIn3[k-1])
                    LcsDynArray[i][j][k] = 1 + LcsDynArray[i-1][j-1][k-1];
                else
                    LcsDynArray[i][j][k] = max(LcsDynArray[i-1][j][k], max(LcsDynArray[i][j-1][k], LcsDynArray[i][j][k-1]));
    cout << LcsDynArray[LcsIn1.size()][LcsIn2.size()][LcsIn3.size()] << endl;
    
    //print the subsequence
    int i = LcsIn1.size(), j = LcsIn2.size(), k = LcsIn3.size();
    vector<char> Lcs3stringOp;
    while(i > 0 && j > 0 && k > 0)
    {
        if(LcsIn1[i-1] == LcsIn2[j-1] && LcsIn2[j-1] == LcsIn3[k-1])
        {
            Lcs3stringOp.push_back(LcsIn1[i-1]);
            --i;
            --j;
            --k;
        }
        else
        {
            if(LcsDynArray[i][j][k] == LcsDynArray[i-1][j][k])
                --i;
            else if(LcsDynArray[i][j][k] == LcsDynArray[i][j-1][k])
                --j;
            else
                --k;
        }
    }
    for(const auto &i : Lcs3stringOp)
        cout << i;
    cout << endl;
    */
    
    /*
    //Longest common substring
    const string LcsubIn1 = "cardiologist", LcsubIn2 = "epidemiologist";
    vector<vector<int>> LcsubDynArray (LcsubIn1.size()+1, vector<int>(LcsubIn2.size()+1,0));
    int maxcnt = INT_MIN, maxi = 0, maxj = 0;
    for(int i = 1; i <= LcsubIn1.size(); ++i)
        for(int j = 1; j <= LcsubIn2.size(); ++j)
            if(LcsubIn1[i-1] == LcsubIn2[j-1])
            {
                LcsubDynArray[i][j] = 1 + LcsubDynArray[i-1][j-1];
                if(maxcnt < LcsubDynArray[i][j])
                {
                    maxcnt = LcsubDynArray[i][j];
                    maxi = i;
                    maxj = j;
                }
            }
    
    cout << maxcnt << " " << maxi << " " << maxj << endl;
    
    //print the subsequence
    int i = maxi, j = maxj;
    vector<char> Lcsub2stringOp;
    while(LcsubDynArray[i][j])
    {
        Lcsub2stringOp.push_back(LcsubIn1[i-1]);
        --i;
        --j;
    }
    for(const auto &i : Lcsub2stringOp)
        cout << i;
    cout << endl;
    */
    
    /*
    //DCp210 - A Collatz sequence in mathematics can be defined as follows. Starting with any positive integer:
    //if n is even, the next number in the sequence is n / 2
      //  if n is odd, the next number in the sequence is 3n + 1
        //    It is conjectured that every such sequence eventually reaches the number 1. Test this conjecture.
          //   Bonus: What input n <= 1000000 gives the longest sequence?
    
    unordered_map<int,int> CollatzMap;
    const int N = 100;
    int MaxCnt = INT_MIN, MaxN = INT_MIN, tempcnt;
    for(int i = 2; i <= N; ++i)
    {
        tempcnt = FindCollatzCnt(i,CollatzMap);
        if(MaxCnt < tempcnt)
        {
            MaxCnt = tempcnt;
            MaxN = i;
        }
    }
    cout << MaxN << " " << MaxCnt << endl;
    */
    
    /*
    //shortest range in 'k' sorted lists
    //Find smallest range containing elements from k lists
    //Given k sorted lists of integers of size n each, find the smallest range that includes at least element from each of the k lists. If more than one smallest ranges are found, print any one of them.
    const vector<vector<int>> Ksortedlist = {{1,20,30,40,80,100,945},{81,85,89,101},{4,45,102,104}};
    vector<int> list_idx(Ksortedlist.size(), 0);
    
    //priority_queue<INT_PAIR,vector<INT_PAIR>,greater<INT_PAIR>> PQ_min_sortlist;
    //priority_queue<INT_PAIR,vector<INT_PAIR>,less<INT_PAIR>> PQ_max_sortlist;
    multiset<INT_PAIR, CompareSortedKlists> MS_sortlist;
    
    //initialize set (binart search tree) with min. elements from the sorted lists
    for(int i = 0; i < Ksortedlist.size(); ++i)
        MS_sortlist.insert(make_pair(Ksortedlist[i][0],i));
    
    int rangemin = INT_MIN, rangemax = INT_MIN, rangediff = INT_MAX, temp_diff;
    INT_PAIR min_elem_idx_pair;
    while(1)
    {
        temp_diff = MS_sortlist.rbegin()->first - MS_sortlist.begin()->first;
        if(temp_diff < rangediff)
        {
            rangediff = temp_diff;
            rangemin = MS_sortlist.begin()->first;
            rangemax = MS_sortlist.rbegin()->first;
        }
        
        //pop the min. element and advance to the next element in the min array's list
        min_elem_idx_pair = *MS_sortlist.begin();
        MS_sortlist.erase(MS_sortlist.begin());
        ++list_idx[min_elem_idx_pair.second];
        if( list_idx[min_elem_idx_pair.second] >= Ksortedlist[min_elem_idx_pair.second].size() )
            break;
        MS_sortlist.insert(make_pair(Ksortedlist[min_elem_idx_pair.second][list_idx[min_elem_idx_pair.second]],min_elem_idx_pair.second));
    }
    cout << rangemin << " " << rangemax << " " << rangediff << endl;
    */
    
    /*
    //DCP212 - Spreadsheets often use this alphabetical encoding for its columns: "A", "B", "C", ..., "AA", "AB", ..., "ZZ", "AAA", "AAB", ....
    //Given a column number, return its alphabetical column id. For example, given 1, return "A". Given 27, return "AA".
    string ExcelCol = "AAZ";
    const int base = 26;
    int ColNum = 0, NumDigits = 0;
    vector<int> DigitVal;
    
    for_each(ExcelCol.begin(), ExcelCol.end(), [](char &i){i=tolower(i);});
    transform(ExcelCol.begin(), ExcelCol.end(), back_inserter(DigitVal),
              [](const char &i){return (static_cast<int>(i)-static_cast<int>('a')+1);});
    
    for(auto i = DigitVal.crbegin(); i != DigitVal.crend(); ++i)
    {
        ColNum += (pow(base,NumDigits)*(*i));
        ++NumDigits;
    }
    cout << ColNum << endl;
    */
    
    /*
    //DCP213 - Given a string of digits, generate all possible valid IP address combinations.
    //IP addresses must follow the format A.B.C.D, where A, B, C, and D are numbers between 0 and 255. Zero-prefixed numbers, such as 01 and 065, are not allowed, except for 0 itself.
    //For example, given "2542540123", you should return ['254.25.40.123', '254.254.0.123'].
    const string Ip = "256000";//2542540123";
    vector<int> TempIp;
    FindIpAddr(Ip, Ip.size()-1, 0, TempIp, 0);
    */
    
    /*
    //DCP214 - Given an integer n, return the length of the longest consecutive run of 1s in its binary representation.
    //For example, given 156, you should return 3.
    int BinNum = 253;
    int OnesCnt = 0, MaxOnesCnt = INT_MIN, BinDigit = 0;
    while(BinNum)
    {
        BinDigit = BinNum & 1;
        BinNum >>= 1;
        if(BinDigit)
            ++OnesCnt;
        else
        {
            if(MaxOnesCnt < OnesCnt)
                MaxOnesCnt = OnesCnt;
            OnesCnt = 0;
        }
    }
    if(MaxOnesCnt < OnesCnt)
        MaxOnesCnt = OnesCnt;
    cout << MaxOnesCnt << endl;
     */
    
    /*
    //DCP216 - Given a number in Roman numeral format, convert it to decimal.
    const string RomanIn = "XLV";
    const unordered_map<char,int> RomanConv = {make_pair('M',1000),make_pair('D',500),make_pair('C',100),
        make_pair('L',50),make_pair('X',10),make_pair('V',5),make_pair('I',1)};
    int DecOutput = 0, PrevVal = INT_MAX, CurrVal = 0;
    auto itr = RomanConv.begin();
    for(const auto &i : RomanIn)
    {
        itr = RomanConv.find(i);
        if(RomanConv.end() == itr)
        {
            cout << "Invalid input " << i << endl;
            DecOutput = INT_MIN;
            break;
        }
        CurrVal = (*itr).second;
        if(CurrVal <= PrevVal)
            DecOutput += CurrVal;
        else
            DecOutput += (CurrVal-PrevVal-PrevVal);
        PrevVal = CurrVal;
    }
    cout << DecOutput << endl;
    */
    
    /*
    //Distance between 2 nodes in a BST
        //Check if the 2 nodes are present in BST
        //Find LCA of the 2 nodes
        //Find distance between each node and LCA, and add them
    const vector<int> TreeInputInt = {4,2,1,3,6,5,7};
    for(const auto &i : TreeInputInt)
        InsertKeyToBst(i, " ");
    PrintBst();
    const int NodeAval = 5, NodeBval = 4;
    shared_ptr<BstType> NodeAPtr = IsNodePresentBst(BstHead, NodeAval);
    shared_ptr<BstType> NodeBPtr = IsNodePresentBst(BstHead, NodeBval);
    if(nullptr != NodeAPtr && nullptr != NodeBPtr)
    {
        shared_ptr<BstType> LcaNodePtr = LcaBst(BstHead, NodeAPtr->key, NodeBPtr->key);
        cout << "BT BST " << (LcaNodePtr ? LcaNodePtr->key : -1) << endl;
        
        cout << BstLevelCnt(LcaNodePtr, NodeAPtr->key) + BstLevelCnt(LcaNodePtr, NodeBPtr->key) << endl;
    }
    */
    
    /*
    //DCP217 - We say a number is sparse if there are no adjacent ones in its binary representation. For example, 21 (10101) is sparse, but 22 (10110) is not. For a given input N, find the smallest sparse number greater than or equal to N.
    //Do this in faster than O(N log N) time.
    const int SparseIn = 38; //1011
    int SparseOut = SparseIn, ContOnesLoc, mask;
    ContOnesLoc = FindContOnesLocFromLsb(SparseOut);
    while(INT_MIN != ContOnesLoc)
    {
        //set everything below the ContOnesLoc to 0, and then set the next bit to 1
        //set all the bits from ContOnesLoc to LSB to 0
        mask = 1 << (ContOnesLoc+1);
        --mask;
        mask = ~mask;
        SparseOut &= mask;
        
        //set (ContOnesLoc + 1) bit to 1
        mask = 1 << (ContOnesLoc+1);
        SparseOut |= mask;
        ContOnesLoc = FindContOnesLocFromLsb(SparseOut);
    }
    cout << SparseOut << endl;
     */
    
    /*
    //compute fibanocci in O(log n) time
    //the top right element (or bottom left element) of matrix ^ N gives the Fib(N) number
    //the matrix is a 2x2 matrix [1,1]
    //                           [1,0]
    const int FibN = 25;
    const vector<vector<int>> Fib2x2array = {{1,1},{1,0}};
    vector<vector<int>> FibOutput;
    ComputeMatrixPower(Fib2x2array, FibN, FibOutput);
    cout << FibOutput[0][1] << endl;
    */

    /*
    //Trie
    shared_ptr<TrieType> TrieHead = make_shared<TrieType>();
    TrieHead->AddWord("banana");
    TrieHead->AddWord("bang");
    TrieHead->AddWord("a");
    TrieHead->AddWord("and");
    TrieHead->AddWord("z");
    TrieHead->AddWord("zoo");
    TrieHead->PrintAllWords();
    */
    
    /*
    //DCP221 - Let's define a "sevenish" number to be one which is either a power of 7, or the sum of unique powers of 7. The first few sevenish numbers are 1, 7, 8, 49, and so on. Create an algorithm to find the nth sevenish number.
    bool ExitLoop = false;
    int pow7;
    vector<int> OpVec;
    auto size = OpVec.size();
    const int SevenN = 20;
    for(int i = 0; i < SevenN; ++i)
    {
        pow7 = static_cast<int>(pow(7, i));
        OpVec.emplace_back(pow7);
        size = OpVec.size() - 1;
        for(int j = 0; j < size; ++j)
        {
            OpVec.emplace_back(OpVec.at(j) + pow7);
            if(SevenN == OpVec.size())
            {
                ExitLoop = true;
                break;
            }
        }
        if(ExitLoop)
            break;
    }
    cout << OpVec.at(SevenN-1) << endl;
    */
    
    /*
    //Multiply 2 integers
    const int num1 = 8826, num2 = 1234;
    int smaller_num = (num1 < num2) ? num1 : num2;
    int larger_num = (num1 < num2) ? num2 : num1;
    
    //Extract each digit of smaller_num without using division/modulo
    ostringstream str_stream;
    str_stream << smaller_num;
    string smaller_num_str = str_stream.str();
    
    int MultAnswer = 0;
    for(int i = 0; i < smaller_num_str.size(); ++i)
    {
        int digit = static_cast<int>(smaller_num_str[smaller_num_str.size()-1-i]) - static_cast<int>('0');
        
        //mulitply larger number by smaller number's digit (by adding larger number 'digit' times)
        int temp = 0;
        for(int j = 0; j < digit; ++j)
            temp += larger_num;
        
        //multiply by power of 10
        MultAnswer += MultByPowerOf10(temp, i);
    }
    cout << MultAnswer << endl;
    */
    
    /*
    //DCP224 - Given a sorted array, find the smallest positive integer that is not the sum of a subset of the array.
    //For example, for the input [1, 2, 3, 10], you should return 7.
    //Do this in O(N) time.
    const vector<int> SortedArray = {1, 2, 3, 4, 6, 10};
    int SubArraySum = SortedArray[0];
    if(SubArraySum > 1)
        cout << "1" << endl;
    else
    {
        for(int i = 1; i < SortedArray.size(); ++i)
        {
            if(SortedArray[i] > SortedArray[i-1] + 1)
                break;
            SubArraySum += SortedArray[i];
        }
        cout << SubArraySum+1 << endl;
    }
    */
    
    /*
    //DCP225 - There are N prisoners standing in a circle, waiting to be executed. The executions are carried out starting with the kth person, and removing every successive kth person going clockwise until there is no one left.
    //Given N and k, write an algorithm to determine where a prisoner should stand in order to be the last survivor.
    //For example, if N = 5 and k = 2, the order of executions would be [2, 4, 1, 5, 3], so you should return 3.
    //Bonus: Find an O(log N) solution if k = 2.
    const int N = 5, K = 2;
    cout << FindKthElemToDelete(N,K) << endl;
    */
    
    /*
    //Word break problem
    //Given a string, break/split into different words available in dictionary
    const string inString = "iamace";
    const unordered_set<string> WordDict = {"i", "am", "ace", "a"};
    cout << SplitSentance(inString, WordDict, 0) << endl;
    */
    
    /*
    //Facotorize a number into primes
    vector<int> PrimeNums = {2}; //look-up table for all prime numbers between 2 and FactNum
    const int FactNum = 439;
    int IsPrime = true;
    for(unsigned int i = 3; i <= FactNum; i=i+2)
    {
        IsPrime = true;
        int SqRoot = sqrt(i);
        for(const auto &j : PrimeNums)
        {
            if(0 == i%j && j <= SqRoot)
            {
                IsPrime = false;
                break;
            }
        }
        if(IsPrime)
            PrimeNums.push_back(i);
    }
    
    vector<int> FactorizeOp;
    const unordered_set<int> PrimeSet (PrimeNums.begin(), PrimeNums.end());
    int idx = 0, temp = FactNum;
    while(temp > 1)
    {
        if(temp % PrimeNums[idx] == 0)
        {
            FactorizeOp.push_back(PrimeNums[idx]);
            temp = temp/PrimeNums[idx];
            if(PrimeSet.find(temp) != PrimeSet.end())
            {
                FactorizeOp.push_back(temp);
                break;
            }
        }
        else
            ++idx;
    }
    for(const auto &i:FactorizeOp)
        cout << i << " ";
    cout << endl;
    */
    
    /*
    //Leetcode 149 - max points on a line
    struct Point
    {
        int x;
        int y;
        Point() : x(0), y(0) {}
        Point(int a, int b) : x(a), y(b) {}
    };
    class Solution
    {
     #if 0
        class LineEq
        {
        public:
            long long int slope;
            long long int constant;
            int num_points;
            Point p1;
            Point p2;
            LineEq(const Point p1, const Point p2) : p1(p1), p2(p2), num_points(2)
            {
                if(p2.x != p1.x)
                {
                    slope = (p2.y-p1.y)/(p2.x-p1.x);
                    constant = p1.y - (slope*p1.x);
                }
                else
                {
                    slope = INT_MAX;
                    constant = p1.x;
                }
            }
            bool IsInLine(const Point p) const
            {
                if(INT_MAX != slope)
                    return ((slope*p.x) + constant == p.y);
                else
                    return (constant == p.x);
            }
        };
     #endif
    public:
        int maxPoints(vector<Point>& points)
        {
            if(points.size() <= 2)
                return points.size();
            
            int MaxCnt = INT_MIN;
            
            for(int i = 0; i < points.size(); ++i)
            {
                unordered_map<long double,int> SlopeSet;
                int Samex = 1;
                int Samepnt = 0;
                for(int j = 0; j < points.size(); ++j)
                {
                    if(i == j)
                        continue;
                    if(points[i].x == points[j].x &&
                       points[i].y == points[j].y)
                    {
                        ++Samepnt;
                    }
                    if(points[i].x == points[j].x)
                    {
                        ++Samex;
                        continue;
                    }
                    long double slope = static_cast<long double>(points[j].y - points[i].y)/static_cast<long double>(points[j].x - points[i].x);
                    if(SlopeSet.find(slope) != SlopeSet.end())
                    {
                        ++SlopeSet[slope];
                        MaxCnt = max(MaxCnt,SlopeSet[slope]);
                    }
                    else
                    {
                        SlopeSet.insert(make_pair(slope, 2+Samepnt));
                        MaxCnt = max(MaxCnt,SlopeSet[slope]);
                    }
                }
                MaxCnt = max(MaxCnt,Samex);
            }
            
            return MaxCnt;
        }
    };
    
    vector<Point> PointsInput;
    PointsInput.emplace_back(0,0);
    PointsInput.emplace_back(94911151,94911150);
    PointsInput.emplace_back(94911152,94911151);
    Solution PointsOnLine;
    cout << PointsOnLine.maxPoints(PointsInput) << endl;
    */

    /*
    //Simple DFS based on array
    //Array is Directed Acyclic Graph (DAG). Every even idx is connected to the next even idx.
    //Every 'multiple of 3' idx is connected to the next 'multiple of 3' idx
    //Perform Topological sort
    vector<char> CharArray;
    const int N = 13;
    for(int i = 0; i < N; ++i)
        CharArray.push_back('A'+i);
    stack<char> TopoSortOp;
    vector<bool> DfsArrayVisited(CharArray.size(), false);
    for(int i = 0; i < N; ++i)
        if(!DfsArrayVisited[i])
            RunDfsArray(i,DfsArrayVisited,TopoSortOp,CharArray);
    while(!TopoSortOp.empty())
    {
        cout << TopoSortOp.top() << " ";
        TopoSortOp.pop();
    }
    cout << endl;
    */
    
    /*
    //Leetcode 753 - cracking the safe
    //There is a box protected by a password. The password is n digits, where each letter can be one of the first k digits 0, 1, ..., k-1.
    //You can keep inputting the password, the password will automatically be matched against the last n digits entered.
    //For example, assuming the password is "345", I can open it when I type "012345", but I enter a total of 6 digits.
    //Please return any string of minimum length that is guaranteed to open the box after the entire string is inputted.
    const int K = 3;
    const int N = 3;
    unordered_set<string> Visited;
    string InitPwd (N,'0');
    Visited.insert(InitPwd);
    string OpPwd = InitPwd;
    cout << CrackSafeDfs(K,N,Visited,OpPwd) << endl;
    cout << OpPwd << endl;
    */
    
    /*
    //DCP228
    //Given a list of numbers, create an algorithm that arranges them in order to form the largest possible integer. For example, given [10, 7, 76, 415], you should return 77641510.
    vector<int> InputVec = {0,0,0,0};
    sort(InputVec.begin(), InputVec.end(), MyCompSortForDcp228);
    
    long long int CombineOp = 0;
    int numDigits = 0, numZero = 0;
    for(auto i = InputVec.crbegin(); i != InputVec.crend(); ++i)
    {
        if(*i == 0)
            ++numZero;
        CombineOp += (*i * pow(10,numDigits));
        numDigits += GetNumberOfDigits(*i);
    }
    CombineOp *= pow(10,numZero);
    cout << CombineOp << endl;
    */
    
    /*
    //DCP229 - Snakes and Ladders is a game played on a 10 x 10 board, the goal of which is get from square 1 to square 100. On each turn players will roll a six-sided die and move forward a number of spaces equal to the result. If they land on a square that represents a snake or ladder, they will be transported ahead or behind, respectively, to a new square.
    //Find the smallest number of turns it takes to play snakes and ladders.
    //For convenience, here are the squares representing snakes and ladders, and their outcomes:
    const unordered_map<int,int> snakes = {make_pair(16,6), make_pair(48,26), make_pair(49,11), make_pair(56,53), make_pair(62,19),
        make_pair(64,60), make_pair(87,24), make_pair(93,73), make_pair(95,75), make_pair(98,78)};
    const unordered_map<int,int> ladders = {make_pair(1,38), make_pair(4,14), make_pair(9,31), make_pair(21,42), make_pair(28,84),
        make_pair(36,44), make_pair(51,67), make_pair(71,91), make_pair(80,100)};
    const int SNL_SIZE = 100;
    vector<bool> snlVisited (SNL_SIZE,false);
    queue<pair<int,int>> snlQ; //1st int indicates the current square position, 2nd int indicates the number of die rolls to get there
    snlQ.push(make_pair(0, 0));
    cout << "push Q 0 0" << endl;
    
    while(!snlQ.empty())
    {
        auto curr = snlQ.front();
        snlQ.pop();
        snlVisited[curr.first] = true;
        cout << "pop Q " << curr.first << " " << curr.second << endl;
        
        if(curr.first == SNL_SIZE-1)
        {
            cout << "Game over in " << curr.second << " steps!!!" << endl;
            break;
        }
        
        //if curr square is a snake, goto the new square without increasing the 'die roll' cnt
        auto isSnake = snakes.find(curr.first);
        if(isSnake != snakes.end())
        {
            pair<int,int> next = make_pair(isSnake->second, curr.second);
            if(!snlVisited[next.first])
            {
                snlVisited[next.first] = true;
                cout << "push Q " << next.first << " " << next.second << endl;
                snlQ.push(next);
            }
            continue;
        }
        
        //if curr square is a ladder, goto the new square without increasing the 'die roll' cnt
        auto isLadder = ladders.find(curr.first);
        if(isLadder != ladders.end())
        {
            pair<int,int> next = make_pair(isLadder->second, curr.second);
            if(!snlVisited[next.first])
            {
                snlVisited[next.first] = true;
                cout << "push Q " << next.first << " " << next.second << endl;
                snlQ.push(next);
            }
            continue;
        }
    
        //at this point, current square neither has a snake or a ladder. Increment the die roll cnt and add to Queue
        for(int i = 1; i <= 6; ++i)
        {
            pair<int,int> next = make_pair(curr.first + i, curr.second + 1);
            if(next.first < SNL_SIZE && !snlVisited[next.first])
            {
                snlVisited[next.first] = true;
                cout << "push Q " << next.first << " " << next.second << endl;
                snlQ.push(next);
            }
        }
    }
    */
    
    /*
    //test code for copy constructor
    TempCl test_class = CreateClass();
    */
    
    /*
    //DCp230 - You are given N identical eggs and access to a building with k floors. Your task is to find the lowest floor that will cause an egg to break, if dropped from that floor. Once an egg breaks, it cannot be dropped again. If an egg breaks when dropped from the xth floor, you can assume it will also break when dropped from any floor greater than x.
    //Write an algorithm that finds the minimum number of trial drops it will take, in the worst case, to identify this floor.
    //For example, if N = 1 and k = 5, we will need to try dropping the egg at every floor, beginning with the first, until we reach the fifth floor, so our solution will be 5.
    const int NUM_EGGS = 3;
    const int NUM_FLOORS = 100;
    unordered_map<string,int> MemoEgg;
    cout << CntNumAttempsEggs(NUM_EGGS, NUM_FLOORS, MemoEgg) << endl;
    */
    
    /*
    //Permute all combinations of a number pad (keypad / key pad)
    //Compute All Mnemonics For A Phone Number
    const string InNum = "7149162206";
    const vector<string> NumAlphTable = {"", "", "ABC", "DEF", "GHI", "JKL", "MNO", "PQRS", "TUV", "WXYZ"};
    vector<string> KeypadOp;
    KeypadOp.reserve(pow(4,InNum.size())); //max will be 4 options per digit and number of digits
    string Tempstr;
    Tempstr.reserve(InNum.size());
    
    FindAllKeypadComb(0, InNum, NumAlphTable, KeypadOp, Tempstr);
    
    for(const auto &i:KeypadOp)
        cout << i << endl;
    */
    
    /*
    //dutch national flag (DNF) problem - sort array of 0s, 1s, and 2s
    vector<int> DNF_in = {0,1,0,2,0,2,1,2,0,0,1,2,2,1,1,0,2,0};
    
    //forward pass - chose value '1' as pivot and move all 0s to left of the array, before 1
    //find 1st location which is >= 1
    int placement_idx0 = 0;
    int PIVOT_VAL = 1;
    while(DNF_in[placement_idx0] < PIVOT_VAL && placement_idx0 < DNF_in.size())
        ++placement_idx0;
    
    //find location which is < pivot, and swap with placement index
    int next_idx = placement_idx0+1;
    while(next_idx < DNF_in.size())
    {
        if(DNF_in[next_idx] < PIVOT_VAL)
        {
            int temp = DNF_in[next_idx];
            DNF_in[next_idx] = DNF_in[placement_idx0];
            DNF_in[placement_idx0] = temp;
            ++placement_idx0;
        }
        ++next_idx;
    }
    
    //backward pass - chose value '1' as pivot and move all 2s to right of the array, after 1
    //do this from the end of the array until placement_idx (elements < placement idx) are 0
    int placement_idx2 = DNF_in.size()-1;
    placement_idx0 %= DNF_in.size();
    while(DNF_in[placement_idx2] > PIVOT_VAL && placement_idx2 >= placement_idx0)
        --placement_idx2;
    
    //find location which is > pivot, and swap with placement index
    next_idx = placement_idx2-1;
    while(next_idx >= 0)
    {
        if(DNF_in[next_idx] > PIVOT_VAL)
        {
            int temp = DNF_in[next_idx];
            DNF_in[next_idx] = DNF_in[placement_idx2];
            DNF_in[placement_idx2] = temp;
            --placement_idx2;
        }
        --next_idx;
    }
    
    for(const auto &i:DNF_in)
        cout << i << " ";
    cout << endl;
    */
    
    /*
    //LRU cache - least recently used
    class LRUcache
    {
    private:
        list<pair<int,int>> itemList; //list that stores the key value pair
        unordered_map<int,list<pair<int,int>>::iterator> itemHash; //hash table that has the iterator (address) from the list of each keys
        const int LRU_size;
    public:
        LRUcache(const int size) : LRU_size(size)
        {
            itemHash.reserve(size);
        }
        int get(const int key)
        {
            if(itemHash.empty())
                return INT_MIN;
            
            auto item_itr = itemHash.find(key);
            if(item_itr == itemHash.end())
                return INT_MIN;
            
            //move the item in the list to the head
            //note: important!!!! - iterator is an address to the node. Moving/splicing the list does not change the iterator
            itemList.splice(itemList.begin(), itemList, item_itr->second);
            
            return itemList.front().second;
        }
        void put(const int key, const int value)
        {
            //check if key already exists
            auto item_itr = itemHash.find(key);
            if(item_itr == itemHash.end())
            {
                //check if cache is full. if so, delete the last used item (tail node)
                if(itemHash.size() == LRU_size)
                {
                    //delete the last item
                    auto itr = itemList.back();
                    itemHash.erase(itr.first);
                    itemList.pop_back();
                }
                
                itemList.push_front(make_pair(key, value));
                itemHash.insert(make_pair(key, itemList.begin()));
            }
            else
            {
                //key already exists. move to head of the list and update to new value
                item_itr->second->second = value;
                itemList.splice(itemList.begin(), itemList, item_itr->second);
            }
        }
    };
    
    LRUcache tempCache(3);
    cout << tempCache.get(3) << endl;
    
    tempCache.put(1, 4);
    tempCache.put(2, 5);
    tempCache.put(3, 6);
    cout << tempCache.get(2) << endl;
    tempCache.put(4, 7);
    tempCache.put(2, 8);
    */
    
    /*
    //test lambda functions
    const vector<int> testVec = {1,2,3,4,5};
    vector<int> mul2Vec;
    for_each(testVec.cbegin(), testVec.cend(), [](int v){cout << v << endl;});
    transform(testVec.cbegin(), testVec.cend(), back_inserter(mul2Vec), [](int v){return v<<1;});
    for_each(mul2Vec.cbegin(), mul2Vec.cend(), [](int v){cout << v << endl;});
    */
    
    /*
    //test lower and upped_bound
    vector<int> testVec = {1,1,2,3,3,3,4,5,5,0,8};
    #define MY_SORT_COMP [](const int a, const int b){return a<b;}
    
    //sort and print output
    sort(testVec.begin(), testVec.end(), MY_SORT_COMP);
    for_each(testVec.cbegin(), testVec.cend(), [](int v){cout << v << endl;});
    
    auto itr1 = lower_bound(testVec.cbegin(), testVec.cend(), 3, MY_SORT_COMP);
    cout << "lower idx " << itr1-testVec.cbegin() << endl;
    
    auto itr2 = upper_bound(testVec.cbegin(), testVec.cend(), 3, MY_SORT_COMP);
    cout << "upper idx " << itr2-testVec.cbegin() << endl;
    
    cout << "binary search " << binary_search(testVec.cbegin(), testVec.cend(), 8, MY_SORT_COMP) << endl;
    
    //find total occureneces of an element 'k' in a sorted array
    const int k = 0;
    auto itr3 = lower_bound(testVec.cbegin(), testVec.cend(), k, MY_SORT_COMP);
    auto itr4 = upper_bound(testVec.cbegin(), testVec.cend(), k, MY_SORT_COMP);
    cout << "total occurences " << itr4-itr3 << endl;
    
    //equal range
    auto itr5 = equal_range(testVec.cbegin(), testVec.cend(), -1, MY_SORT_COMP);
    for(auto i = itr5.first; i != itr5.second; ++i)
        cout << *i << endl;
    */
    
    /*
    //DCP231 - Given a string with repeated characters, rearrange the string so that no two adjacent characters are the same. If this is not possible, return None.
    //For example, given "aaabbc", you could return "ababac". Given "aaab", return None.
    const string RepeatChar = "abcc";
    list<char> CharList (RepeatChar.cbegin(), RepeatChar.cend());
    string OpStr;
    OpStr.reserve(RepeatChar.size());
    cout << RearrangeRepeatChar(CharList, OpStr) << " ";
    cout << OpStr << endl;
    */
    
    /*
    //DCP231 - Given a string with repeated characters, rearrange the string so that no two adjacent characters are the same. If this is not possible, return None.
    //For example, given "aaabbc", you could return "ababac". Given "aaab", return None.
    //using pq
    const string RepeatChar = "aaab";
    unordered_map<char,int> CharCnt; //256 ascii characters
    for(const auto & i : RepeatChar)
        ++CharCnt[i];
#define PQ_REARRANGE pair<int,char>
    class PQ_REARRANGE_COMP
    {
    public:
        bool operator()(const PQ_REARRANGE a, const PQ_REARRANGE b)
        {
            return a.first < b.first;
        }
    };
    priority_queue<PQ_REARRANGE, vector<PQ_REARRANGE>, PQ_REARRANGE_COMP> RearrangePQ;
    for(const auto & i : CharCnt)
        RearrangePQ.push(make_pair(i.second, i.first));
    
    string OpStr;
    OpStr.reserve(RepeatChar.size());
    while(!RearrangePQ.empty())
    {
        //greedy algm - push the most frequent element from the PQ and remove it
        auto PQ_top = RearrangePQ.top();
        RearrangePQ.pop();
        
        if(OpStr.empty() || OpStr.back() != PQ_top.second)
        {
            OpStr.push_back(PQ_top.second);
            --PQ_top.first;
        }
        else
        {
            //put the next character and then add it back to Q
            if(RearrangePQ.empty())
                break;
            
            auto PQ_2ndtop = RearrangePQ.top();
            RearrangePQ.pop();
            OpStr.push_back(PQ_2ndtop.second);
            --PQ_2ndtop.first;
            if(PQ_2ndtop.first)
                RearrangePQ.push(PQ_top);
        }
        
        //if the 1st top element is > 0, add it back to Q
        if(PQ_top.first)
            RearrangePQ.push(PQ_top);
    }
    
    if(OpStr.size() == RepeatChar.size())
        cout << "sucess!! " << OpStr << endl;
    else
        cout << "cannot rearrange " << OpStr << endl;
    */
    
    //DCP219 - connect 4 game - Connect 4 is a game where opponents take turns dropping red or black discs into a 7 x 6 vertically suspended grid. The game ends either when one player creates a line of four consecutive discs of their color (horizontally, vertically, or diagonally), or when there are no more spots left in the grid.
    //array 7x6; 0 - empty slot, 1 - red, 2 - yellow
    
    /*
    //Total ways to decode a number/string. For example, 112 can be AAB, KB, AL (A = 1, B = 2, ..., Z = 26)
    const int DecodeNum = 1126;
    int ans = 0;
    FindDecodeNumWays(DecodeNum, ans);
    cout << ans << endl;
    */
    
    /*
    //Structure zero padding
    struct test {
        char a;
        long long int x;
        char b;
    };
    cout << sizeof(test) << endl;
     */
    
    /*
    //Convert 2 integers (int and frac portion) into binary format (IEEE 754)
    //1 bit for sign
    //8 bits for exponent (scaled/biased by 2^127)
    //23 bits of mantissa
    const int int_val = 48;
    const int frac_val = 4501; //read as 0.4501
    vector<int> Int_binary, Frac_binary;
    
    //convert integer to binary - keep dividing by 2 until quotient is 0. For every division, the remainder will
    //be the binary representation
    int temp = int_val;
    while(temp)
    {
        Int_binary.push_back(temp % 2);
        temp /= 2;
    }
    reverse(Int_binary.begin(), Int_binary.end());
    
    //convert fractional portion - keep multiplying by 2 (for N precision bits). Everytime the result is >= 1, subtract 1
    //and continue multiplying by 2 for N times. The integer portion of the multiplied output (0 or 1) will be the binary
    //representation
    temp = frac_val;
    int NumFracDigits = floor(log10(frac_val)) + 1;
    int ThreshGreaterThan1 = pow(10,NumFracDigits);
    int BinaryPrecisionBits = 31;
    for(int i = 0; i < BinaryPrecisionBits; ++i)
    {
        temp <<= 1;
        if(temp >= ThreshGreaterThan1)
        {
            Frac_binary.push_back(1);
            temp -= ThreshGreaterThan1;
        }
        else
            Frac_binary.push_back(0);
    }
    
    for(const auto &i:Int_binary)
        cout << i;
    cout << ".";
    for(const auto &i:Frac_binary)
        cout << i;
    cout << endl;
    
    int sign = 0;
    int exponent = 127 + Int_binary.size()-1; //move the decimal point to the left (multiply by 2^num of bits moved)
    int mantissa; // mantissa will be the remaining bits of Int_binary (except the MSB) + all the bits of Frac_binary and cut to a length of 23 
    */
    
    /*
    //EPI - 4.6 - find quotient without using arithmetic operations
    const int dividend = 8768764, divisor = 1;
    cout << FindQuotient(dividend, divisor) << endl;
    */
    
    /*
    //LIS - longest increasing subsequence - O(nlogn) solution
    const vector<int> LisInput = {0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15};
    vector<int> EndElem (LisInput.size(), INT_MAX);
    int temp_index;
    for(const auto &i:LisInput)
    {
        //find the earliest index where the new element 'i' can be placed (this is the end element)
        temp_index = lower_bound(EndElem.cbegin(), EndElem.cend(), i) - EndElem.cbegin();
        EndElem[temp_index] = i;
    }
    cout << lower_bound(EndElem.cbegin(), EndElem.cend(), INT_MAX) - EndElem.cbegin() << endl;
    
    //print the elements
    for(const auto &i:EndElem)
        if(i != INT_MAX)
            cout << i << " ";
        else
            break;
    cout << endl;
    */
    
    /*
    //DCP241 - In academia, the h-index is a metric used to calculate the impact of a researcher's papers. It is calculated as follows:
    //A researcher has index h if at least h of her N papers have h citations each. If there are multiple h satisfying this formula, the maximum is chosen.
    //For example, suppose N = 5, and the respective citations of each paper are [4, 3, 0, 1, 5]. Then the h-index would be 3, since the researcher has 3 papers with at least 3 citations.
    //Given a list of paper citations of a researcher, calculate their h-index.
    vector<int> HidxInput = {1,4,1,4,2,1,3,5,6};//{4, 3, 0, 1, 5};
    sort(HidxInput.begin(), HidxInput.end());
    int hidx = 0;
    for(int i = HidxInput.size()-1; i >= 0; --i)
    {
        if(HidxInput.size() - i >= HidxInput[i])
        {
            hidx = HidxInput.size() - i;
            break;
        }
    }
    cout << hidx << endl;
    */
    
    /*
    //DCP245 - You are given an array of integers, where each element represents the maximum number of steps that can be jumped going forward from that element. Write a function to return the minimum number of jumps you must take in order to get from the start to the end of the array.
    //For example, given [6, 2, 4, 0, 5, 1, 1, 4, 2, 9], you should return 2, as the optimal solution involves jumping from 6 to 5, and then from 5 to 9.
    const vector<int> StepsIn = {6, 2, 4, 0, 1, 1, 1, 4, 2, 9};
    vector<int> MaxReach (StepsIn.size(), 0);
    for(int i = 0; i < StepsIn.size(); ++i)
        MaxReach[i] = i + StepsIn[i];
    
    int StepsCounter = 1;
    int temp = EarliestReachableStep(MaxReach, 0, MaxReach.size()-1, MaxReach.size()-1);
    while((temp != 0) && (temp != INT_MAX))
    {
        ++StepsCounter;
        cout << temp << " ";
        temp = EarliestReachableStep(MaxReach, 0, temp, temp);
    }
    if(temp == INT_MAX)
        cout << "\ncannot reach end" << endl;
    else
        cout << "\n" << StepsCounter << endl;
    */
    
    /*
    //DCP246 - Given a list of words, determine whether the words can be chained to form a circle. A word X can be placed in front of another word Y in a circle if the last character of X is same as the first character of Y.
    //For example, the words ['chair', 'height', 'racket', touch', 'tunic'] can form the following circle: chair --> racket --> touch --> height --> tunic --> chair.
    const vector<string> ChainWords = {"chair", "height", "racket", "touch", "tunic"};
    stack<pair<string,int>> WordStack;
    vector<bool> Visited (ChainWords.size(), false);
    vector<int> Indegree (ChainWords.size(), 0);
    vector<int> Outdegree (ChainWords.size(), 0);
    
    WordStack.push(make_pair(ChainWords[0], 0));
    string curr_word;
    int curr_idx;
    
    while(!WordStack.empty())
    {
        curr_word = WordStack.top().first;
        curr_idx = WordStack.top().second;
        WordStack.pop();
        Visited[curr_idx] = true;
        
        for(int i = 0; i < ChainWords.size(); ++i)
        {
            if( (i != curr_idx) && (curr_word.back() == ChainWords[i].front()) )
            {
                ++Outdegree[curr_idx];
                ++Indegree[i];
                if(!Visited[i])
                {
                    WordStack.push(make_pair(ChainWords[i], i));
                }
            }
        }
    }
    */
    
    /*
    //DCP248 - Find the maximum of two numbers without using any if-else statements, branching, or direct comparisons.
    const int x=-100, y=-100;
    cout << (x^(x^y)&((x-y)>>31)) << endl;
    */
    
    /*
    //test smart pointer destructor
    class A;
    class A
    {
        int x;
        shared_ptr<A> next;
    public:
        A() : x(0), next(nullptr)
        {
            cout << "default ctor " << endl;
        }
        A(const int y) : x(y), next(nullptr)
        {
            cout << "value ctor " << x << endl;
        }
        ~A()
        {
            //next = nullptr;
            cout << "dtor " << x << endl;
        }
        shared_ptr<A>& GetNext(void)
        {
            return next;
        }
    };
    shared_ptr<A> Aptr = make_shared<A>();//, Bptr = make_shared<A>(2);
    Aptr->GetNext() = make_shared<A>(2);
    Aptr = nullptr;
    cout << "terminating program" << endl;
    */
    
    /*
    //Given two strings - S1 and S2.
    //Arrange the characters of S1 in same alphabetical order as the characters of S2.
    //If a character of S1 is not present in S2 - such characters should come at the end of the result string, but make sure to retain the order of such characters
    //Case sensitivity is irrelevant
    //e.g. S1 = "Google", S2 = "dog"
    //Output = "ooggle"
    
    //e.g. S1 = "abcdedadf", S2 = "cae"
    //Output = "caaebdddf"
    string S1 = "abcdedadf";
    sort(S1.begin(), S1.end(), MyStringSortComp);
    cout << S1 << endl;
    */
    
    /*
    //Radix sort - similar to counting sort (or LSD sort for strings)
    vector<int> In = {45,90,34,135,2,0,1,45,56,29,21,39,93,2980,1043};
    const int BaseDigit = 10;
    int MaxNum = INT_MIN;
    for(const auto &i : In)
        MaxNum = max(MaxNum,i);
    const int MaxDigits = log10(MaxNum)+1;
    vector<int> AuxArray (In.size(),0);
    for(int i = 0; i < MaxDigits; ++i)
    {
        vector<int> Cnt (BaseDigit+1,0);
        int currdigit;
        
        //cnt the number of numbers with the same digit
        for(int j = 0; j < In.size(); ++j)
        {
            currdigit = static_cast<int>(In[j]/pow(BaseDigit,i)) % BaseDigit;
            ++Cnt[currdigit+1];
        }
        
        //update the cnt array to get the array index for aux vector
        for(int j = 1; j < Cnt.size(); ++j)
            Cnt[j] += Cnt[j-1];
        
        //copy data into aux array based
        for(int j = 0; j < In.size(); ++j)
        {
            currdigit = static_cast<int>(In[j]/pow(BaseDigit,i)) % BaseDigit;
            AuxArray[ Cnt[currdigit]++ ] = In[j];
        }
        
        //copy data back into input array
        for(int j = 0; j < AuxArray.size(); ++j)
            In[j] = AuxArray[j];
    }
    
    //print sorted data
    for(const auto &i : In)
        cout << i << " ";
    cout << endl;
    */
    
    /*
    //generate gray code - can be used to implement towers of hanoi
    const int Nbits = 3;
    vector<int> graycode (pow(2,Nbits),0);
    for(int i = 1; i <= Nbits; ++i)
    {
        int offset = pow(2,i-1);
        for(int j = 0; j < offset; ++j)
            graycode[offset+j] = offset + graycode[offset-j-1];
    }
    
    for(const auto &i:graycode)
        cout << bitset<Nbits>(i) << endl;
    */
    
    /*
    //leetcode 79 - Given a 2D board and a word, find if the word exists in the grid.
    //The word can be constructed from letters of sequentially adjacent cell, where "adjacent" cells are those horizontally or vertically neighboring. The same letter cell may not be used more than once.
    //board =
    //[
     //['A','B','C','E'],
     //['S','F','C','S'],
     //['A','D','E','E']
     //]
    //Given word = "ABCCED", return true.
    //Given word = "SEE", return true.
    //Given word = "ABCB", return false.
    class Solution {
    public:
        bool exist(vector<vector<char>>& board, string word) {
            bool found = false;
            vector<vector<bool>> visited (board.size(), vector<bool>(board[0].size(), false));
            
            for(int i = 0; i < board.size(); ++i)
            {
                for(int j = 0; j < board[i].size(); ++j)
                {
                    if(WordSearchDfs(board, i, j, word, 0, visited))
                    {
                        found = true;
                        break;
                    }
                }
                if(found)
                    break;
            }
            return found;
        }
        
        bool WordSearchDfs(const vector<vector<char>>& board, const int row, const int col, const string &word, const int word_idx, vector<vector<bool>> &visited)
        {
            if(word_idx >= word.size())
                return true;
            
            if(board[row][col] == word[word_idx])
            {
                visited[row][col] = true;
                if(row+1 < board.size() && !visited[row+1][col])
                {
                    if(WordSearchDfs(board, row+1, col, word, word_idx+1, visited))
                        return true;
                }
                if(row-1 >= 0 && !visited[row-1][col])
                {
                    if(WordSearchDfs(board, row-1, col, word, word_idx+1, visited))
                        return true;
                }
                if(col+1 < board[row].size() && !visited[row][col+1])
                {
                    if(WordSearchDfs(board, row, col+1, word, word_idx+1, visited))
                        return true;
                }
                if(col-1 >= 0 && !visited[row][col-1])
                {
                    if(WordSearchDfs(board, row, col-1, word, word_idx+1, visited))
                        return true;
                }
                
                if(word_idx+1 >= word.size())
                    return true;
            }
            
            visited[row][col] = false;
            return false;
        }
    };
    */
    
    /*
    //DCP257 - Given an array of integers out of order, determine the bounds of the smallest window that must be sorted in order for the entire array to be sorted. For example, given [3, 7, 5, 6, 9], you should return (1, 3).
    const vector<int> In = {1,2,3,4,1};
    int leftidx = INT_MIN, rightidx = INT_MIN;
    int leftmax = INT_MIN, rightmin = INT_MAX;
    for(int i = 0; i < In.size(); ++i)
    {
        if(In[i] < leftmax)
            rightidx = i;
        leftmax = max(leftmax, In[i]);
    }
    for(int i = In.size()-1; i >= 0; --i)
    {
        if(In[i] > rightmin)
            leftidx = i;
        rightmin = min(rightmin, In[i]);
    }
    cout << leftidx << " " << rightidx << endl;
    */
    
    /*
    //DCP258 - In Ancient Greece, it was common to write text with the first line going left to right, the second line going right to left, and continuing to go back and forth. This style was called "boustrophedon".
    //Given a binary tree, write an algorithm to print the nodes in boustrophedon order.
    //For example, given the following tree:
    //      1
    //   /     \
    //  2       3
    // / \     / \
    //4   5   6   7
    //You should return [1, 3, 2, 4, 5, 6, 7].
    const vector<int> BtInputInt = {1,2,3,4,6,5,7};
    for(const auto &i : BtInputInt)
        InsertKeyToAnyTree(i, "a");
    PrintBst();
    bool l_to_r = true;
    stack<shared_ptr<BstType>> ltor_s, rtol_s;
    ltor_s.push(BstHead);
    shared_ptr<BstType> curr;
    while(!ltor_s.empty() || !rtol_s.empty())
    {
        if(l_to_r)
        {
            curr = ltor_s.top();
            ltor_s.pop();
            cout << curr->key << " ";
            if(curr->left)
                rtol_s.push(curr->left);
            if(curr->right)
                rtol_s.push(curr->right);
            if(ltor_s.empty())
                l_to_r = false;
        }
        else
        {
            curr = rtol_s.top();
            rtol_s.pop();
            cout << curr->key << " ";
            if(curr->right)
                ltor_s.push(curr->right);
            if(curr->left)
                ltor_s.push(curr->left);
            if(rtol_s.empty())
                l_to_r = true;
        }
    }
    cout << endl;
    */
    
    /*
    //DCP273 - A fixed point in an array is an element whose value is equal to its index. Given a sorted array of distinct elements, return a fixed point, if one exists. Otherwise, return False.
    //For example, given [-6, 0, 2, 40], you should return 2. Given [1, 5, 7, 8], you should return False.
    const vector<int> In = {1, 5, 7, 8};//{-6, 0, 2, 40};
    cout << FixedPntArray(In, 0, In.size()-1) << endl;
    */
    
    /*
    //DCP282 - Given an array of integers, determine whether it contains a Pythagorean triplet. Recall that a Pythogorean triplet (a, b, c) is defined by the equation a2+ b2= c2.
    const vector<int> In = {3,1,4,8,4,9,-9,5,2,-2};
    unordered_multiset<int> SqSet;
    SqSet.reserve(In.size());
    for(auto i = In.cbegin(); i != In.cend(); ++i)
        SqSet.insert((*i) * (*i));
    bool TripFound = false;
    
    for(int i = 0; i < In.size(); ++i)
    {
        auto elem1 = In[i] * In[i];
        auto itr1 = SqSet.equal_range(elem1);
        SqSet.erase(itr1.first);
        for(int j = i+1; j < In.size(); ++j)
        {
            auto elem2 = In[j] * In[j];
            auto itr2 = SqSet.equal_range(elem2);
            SqSet.erase(itr2.first);
            if(SqSet.find(elem1+elem2) != SqSet.end())
            {
                cout << In[i] << " " << In[j] << " " << elem1+elem2 << endl;
                TripFound = true;
                break;
            }
            SqSet.insert(elem2);
        }
        SqSet.insert(elem1);
        if(TripFound)
            break;
    }
    cout << TripFound << endl;
    */
    
    /*
    //DCP283 - A regular number in mathematics is defined as one which evenly divides some power of 60. Equivalently, we can say that a regular number is one whose only prime divisors are 2, 3, and 5.
    //These numbers have had many applications, from helping ancient Babylonians keep time to tuning instruments according to the diatonic scale.
    //Given an integer N, write a program that returns, in order, the first N regular numbers.
    const int N = 20;
    vector<int> RegNumOp;
    RegNumOp.reserve(N);
    int cnt2, cnt3, cnt5;
    cnt2 = cnt3 = cnt5 = 1;
    RegNumOp.push_back(2*cnt2 * 3*cnt3 * 5*cnt5);
    while(RegNumOp.size() < N)
    {
        int mult2 = 2*(cnt2+1) * 3*cnt3 * 5*cnt5;
        int mult3 = 2*cnt2 * 3*(cnt3+1) * 5*cnt5;
        int mult5 = 2*cnt2 * 3*cnt3 * 5*(cnt5+1);
        int minnum = min(mult2, mult3);
        minnum = min(minnum, mult5);
        RegNumOp.push_back(minnum);
        if(minnum == mult2)
            ++cnt2;
        else if(minnum == mult3)
            ++cnt3;
        else
            ++cnt5;
    }
    for(auto i = RegNumOp.cbegin(); i != RegNumOp.cend(); ++i)
        cout << *i << " ";
    cout << endl;
    */
    
    /*
    //DCP290 - On a mysterious island there are creatures known as Quxes which come in three colors: red, green, and blue. One power of the Qux is that if two of them are standing next to each other, they can transform into a single creature of the third color.
    //Given N Quxes standing in a line, determine the smallest number of them remaining after any possible sequence of such transformations.
    //For example, given the input ['R', 'G', 'B', 'G', 'B'], it is possible to end up with a single Qux through the following steps:
    //
    //Arrangement       |   Change
    //----------------------------------------
    //['R', 'G', 'B', 'G', 'B'] | (R, G) -> B
    //['B', 'B', 'G', 'B']      | (B, G) -> R
    //['B', 'R', 'B']           | (R, B) -> G
    //['B', 'G']                | (B, G) -> R
    //['R']                     |
    //R = 1, G = 2, B = 3
    list<int> RgbIn = {1,1,3,2};//{1,2,3,2,3};
    cout << QuxRgb(RgbIn) << endl;
    */
    
    //DCP293 - You have N stones in a row, and would like to create from them a pyramid. This pyramid should be constructed such that the height of each stone increases by one until reaching the tallest stone, after which the heights decrease by one. In addition, the start and end stones of the pyramid should each be one stone high.
    //You can change the height of any stone by paying a cost of 1 unit to lower its height by 1, as many times as necessary. Given this information, determine the lowest cost method to produce this pyramid.
    //For example, given the stones [1, 1, 3, 3, 2, 1], the optimal solution is to pay 2 to create [0, 1, 2, 3, 2, 1].
    
    /*
    //DCP292 - A teacher must divide a class of students into two teams to play dodgeball. Unfortunately, not all the kids get along, and several refuse to be put on the same team as that of their enemies.
    //Given an adjacency list of students and their enemies, write an algorithm that finds a satisfactory pair of teams, or returns False if none exists.
    //For example, given the following enemy graph you should return the teams {0, 1, 4, 5} and {2, 3}.
    //students = {
    //    0: [3],
    //    1: [2],
    //    2: [1, 4],
    //    3: [0, 4, 5],
    //    4: [2, 3],
    //    5: [3] }
    const vector<vector<int>> Enemy = {{3},{2},{1,4},{0,4,5},{2,3},{3}};
    unordered_set<int> team1, team2;
    queue<int> q1,q2;
    vector<bool> visited (Enemy.size(), false);
    
    //put the 1st student in q1
    q1.push(0);
    
    bool currteam1 = true;
    int curr_stud;
    
    while(1)
    {
        if(currteam1)
        {
            curr_stud = q1.front();
            q1.pop();
        
            visited[curr_stud] = true;
            team1.insert(curr_stud);
            for(int i = 0; i < Enemy[curr_stud].size(); ++i)
            {
                if(!visited[Enemy[curr_stud][i]])
                    q2.push(Enemy[curr_stud][i]);
            }
            
            if(q1.empty())
                currteam1 = false;
        }
        else
        {
            curr_stud = q2.front();
            q2.pop();
            
            visited[curr_stud] = true;
            team2.insert(curr_stud);
            for(int i = 0; i < Enemy[curr_stud].size(); ++i)
            {
                if(!visited[Enemy[curr_stud][i]])
                    q1.push(Enemy[curr_stud][i]);
            }
            
            if(q2.empty())
                currteam1 = true;
        }
        
        if(q1.empty() && q2.empty())
            break;
    }
    
    //check if there are any overlaps
    for(auto i = team1.cbegin(); i != team1.end(); ++i)
        if(team2.find( *i ) != team2.end())
        {
            cout << "cannot split into teams!" << *i << endl;
        }
    */
    
    /*
    //DCP295 - Pascal's triangle is a triangular array of integers constructed with the following formula:
    //The first row consists of the number 1.
    //For each subsequent row, each element is the sum of the numbers directly above it, on either side.
    //For example, here are the first few rows:
    //    1
    //   1 1
    //  1 2 1
    // 1 3 3 1
    //1 4 6 4 1
    //Given an input k, return the kth row of Pascal's triangle. Bonus: Can you do this using only O(k) space?
    const int pascal_size = 5;
    vector<int> curr, prev;
    curr.reserve(pascal_size);
    prev.reserve(pascal_size);
    prev.push_back(1);
    for(int i = 1; i < pascal_size; ++i)
    {
        curr.push_back(1);
        for(int j = 1; j < prev.size(); ++j)
            curr.push_back(prev[j] + prev[j-1]);
        curr.push_back(1);
        prev.clear();
        prev = curr;
        curr.clear();
        
        //debug print
        for(auto j = prev.cbegin(); j != prev.cend(); ++j)
            cout << *j << " ";
        cout << endl;
    }
    */
    
    /*
    //DCP297 - At a popular bar, each customer has a set of favorite drinks, and will happily accept any drink among this set. For example, in the following situation, customer 0 will be satisfied with drinks 0, 1, 3, or 6.
    //preferences = {
    //    0: [0, 1, 3, 6],
    //    1: [1, 4, 7],
    //    2: [2, 4, 7, 5],
    //    3: [3, 2, 5],
    //    4: [5, 8]
    //}
    //A lazy bartender working at this bar is trying to reduce his effort by limiting the drink recipes he must memorize. Given a dictionary input such as the one above, return the fewest number of drinks he must learn in order to satisfy all customers.
    //For the input above, the answer would be 2, as drinks 1 and 5 will satisfy everyone.
    const vector<vector<int>> CustIn = {{0,1,3,6},{1,4,7},{2,4,7,5},{3,2,5},{5,8}};
    unordered_set<int> AllDrinks;
    
    for(int i = 0; i < CustIn.size(); ++i)
        for(int j = 0; j < CustIn[i].size(); ++j)
            AllDrinks.insert(CustIn[i][j]);
    
    vector<vector<int>> DrinksList;
    const int NumDrinks = AllDrinks.size();
    for(int i = 0; i < NumDrinks; ++i)
        DrinksList.emplace_back(vector<int>{});
    
    for(int i = 0; i < CustIn.size(); ++i)
        for(int j = 0; j < CustIn[i].size(); ++j)
            DrinksList[ CustIn[i][j] ].push_back(i);
    
    unordered_set<int> AllCust;
    unordered_set<int> BestDrink;
    for(int i = 0; i < CustIn.size(); ++i)
        AllCust.insert(i);
    const int NumCust = CustIn.size();
    
    vector<unordered_set<int>> MinDrinks;
    MinDrinks.emplace_back(set<int>{});
    
    for(int i = 0; i < NumDrinks; ++i)
    {
        CheckCustDrinksComb(i, AllCust, MinDrinks, DrinksList, BestDrink);
        for(int j = 0; j < CustIn.size(); ++j)
            AllCust.insert(i);
    }
    */
    
    /*
    //DCP300 - On election day, a voting machine writes data in the form (voter_id, candidate_id) to a text file. Write a program that reads this file as a stream and returns the top 3 candidates at any given time. If you find a voter voting more than once, report this as fraud.
    ifstream votingmachine("/Users/hari/C-programs/Xcode/TestCPP/TestCPP/TestCPP/vot_machine.txt");
    
    int voterid, candidateid;
    unordered_set<int> voterset;
    unordered_map<int,int> candidate_cnt;
    while((votingmachine >> voterid) && (votingmachine >> candidateid))
    {
        if( !voterset.insert(voterid).second )
        {
            //voter has voted twice (already present in voter list)
            cout << "error. voter already voted" << voterid << endl;
            continue;
        }
        
        candidate_cnt[candidateid]++;
    }
    
    votingmachine.close();
    
    //select top 3 candidates
    vector<pair<int,int>> topcandidate;
    topcandidate.reserve(candidate_cnt.size());
    
    for(auto i = candidate_cnt.cbegin(); i != candidate_cnt.cend(); ++i)
        topcandidate.push_back(make_pair((*i).first, (*i).second));
    
    sort(topcandidate.begin(), topcandidate.end(), MyCandSort);
    */
    
    /*
    //DCP314 - You are the technical director of WSPT radio, serving listeners nationwide. For simplicity's sake we can consider each listener to live along a horizontal line stretching from 0 (west) to 1000 (east).
    //Given a list of N listeners, and a list of M radio towers, each placed at various locations along this line, determine what the minimum broadcast range would have to be in order for each listener's home to be covered.
    //For example, suppose listeners = [1, 5, 11, 20], and towers = [4, 8, 15]. In this case the minimum range would be 5, since that would be required for the tower at position 15 to reach the listener at position 20.
    vector<int> towers = {15,4,8};
    const vector<int> listeners = {20,1,11,5};
    
    //sort list of towers
    sort(towers.begin(), towers.end());
    
    //for all listeners, find the closest tower and keep track of the max. distance
    int maxdist = INT_MIN, mindist = INT_MAX;
    int listidx = 0;
    for(int i = 0; i < listeners.size(); ++i)
    {
        listidx = FindPosSortedArray(towers, 0, towers.size()-1, listeners[i]);
        mindist = INT_MAX;
        
        mindist = min( mindist, abs(listeners[i] - towers[listidx]) );
        if(listidx > 0)
            mindist = min( mindist, abs(listeners[i] - towers[listidx-1]) );
        if(listidx < towers.size()-1)
            mindist = min( mindist, abs(listeners[i] - towers[listidx+1]) );
        
        maxdist = max( maxdist, mindist );
    }
    cout << maxdist << endl;
    */
    
    /*
    //DCP313 - You are given a circular lock with three wheels, each of which display the numbers 0 through 9 in order. Each of these wheels rotate clockwise and counterclockwise.
    //In addition, the lock has a certain number of "dead ends", meaning that if you turn the wheels to one of these combinations, the lock becomes stuck in that state and cannot be opened.
    //Let us consider a "move" to be a rotation of a single wheel by one digit, in either direction. Given a lock initially set to 000, a target combination, and a list of dead ends, write a function that returns the minimum number of moves required to reach the target state, or None if this is impossible.
    vector<int> NextComb = GetLockNextComb(999, 3);
    vector<bool> Visited (1000, false);
    unordered_set<int> Deadends = {554,556,545,565,455,655};
    const int Dest = 555;
    int NumSteps = 0, CurrComb;
    
    //mark all deadends as already visited
    for(auto i = Deadends.cbegin(); i != Deadends.cend(); ++i)
        Visited[*i] = true;
    
    queue<int> q1, q2;
    queue<int> *Currq = &q1;
    bool UsingQ1 = true;
    q1.push(0);
    
    while(!q1.empty() || !q2.empty())
    {
        CurrComb = Currq->front();
        Currq->pop();
        Visited[CurrComb] = true;
        
        if(CurrComb == Dest)
        {
            cout << "Reached Dest!!! " << NumSteps << endl;
            break;
        }
        
        NextComb = GetLockNextComb(CurrComb, 3);
        for(const auto &i : NextComb)
        {
            if(!Visited[i])
            {
                if(UsingQ1)
                    q2.push(i);
                else
                    q1.push(i);
            }
        }
        
        if(Currq->empty())
        {
            if(q1.empty())
            {
                UsingQ1 = false;
                Currq = &q2;
            }
            else
            {
                UsingQ1 = true;
                Currq = &q1;
            }
            ++NumSteps;
        }
    }
    */
    
    return 0;
}

vector<int> GetLockNextComb(const int curr, const int NumDigits)
{
    vector<int> Comb, IntDigits;
    Comb.reserve(NumDigits * 2);
    IntDigits.reserve(NumDigits);
    
    int digit, newnum;
    for(int i = 0; i < NumDigits; ++i)
    {
        digit = curr / static_cast<int>(pow(10,i));
        digit %= 10;
        IntDigits.push_back(digit);
    }
    
    for(int i = 0; i < IntDigits.size(); ++i)
    {
        //reduce digit by 1
        --IntDigits[i];
        if(IntDigits[i] < 0)
            IntDigits[i] += 10;
        
        newnum = 0;
        for(int j = 0; j < IntDigits.size(); ++j)
            newnum += (static_cast<int>(pow(10,j)) * IntDigits[j]);
        Comb.push_back(newnum);
     
        //increase digit by 2
        IntDigits[i] += 2;
        IntDigits[i] %= 10;
        
        newnum = 0;
        for(int j = 0; j < IntDigits.size(); ++j)
            newnum += (static_cast<int>(pow(10,j)) * IntDigits[j]);
        Comb.push_back(newnum);
        
        //goback to original digit (reduce by 1)
        --IntDigits[i];
        if(IntDigits[i] < 0)
            IntDigits[i] += 10;
    }
    
    return Comb;
}

int FindPosSortedArray(const vector<int> &In, const int begin, const int end, const int num)
{
    if(begin >= end)
        return begin;
    
    int mid = begin+end;
    mid >>= 1;
    
    if(num <= In[mid])
        return FindPosSortedArray(In, begin, mid, num);
    return FindPosSortedArray(In, mid+1, end, num);
}

bool MyCandSort(const pair<int,int> &A, const pair<int,int> &B)
{
    return (A.second < B.second);
}

void CheckCustDrinksComb(int currdrink, unordered_set<int> &CustSet, vector<unordered_set<int>> &MinDrinks, const vector<vector<int>> &DrinksList, unordered_set<int> &BestDrink)
{
    unordered_set<int> RemovedCustList;
    
    //remove list of all customers that have this drink
    for(int i = 0; i < DrinksList[currdrink].size(); ++i)
        if( CustSet.erase(DrinksList[currdrink][i]) )
        {
            RemovedCustList.insert(DrinksList[currdrink][i]);
            
            //atleast one more customer is removed. add current drink to list and proceed to next combination
            BestDrink.insert(currdrink);
        }
    
    if(CustSet.empty())
    {
        MinDrinks.emplace_back(BestDrink);
        
        //add list of all customers that have this drink back to the set
        for(auto i = RemovedCustList.cbegin(); i != RemovedCustList.cend(); ++i)
            CustSet.insert(*i);
        
        //Remove this drink from list and proceed to next combination
        BestDrink.erase(currdrink);
        
        return;
    }
    
    //Try the remaining list of drinks
    for(int i = currdrink+1; i < DrinksList.size(); ++i)
    {
        CheckCustDrinksComb(i, CustSet, MinDrinks, DrinksList, BestDrink);
    }
    
    //add list of all customers that have this drink back to the set
    for(auto i = RemovedCustList.cbegin(); i != RemovedCustList.cend(); ++i)
    {
        CustSet.insert(*i);
    
        //Remove this drink from list and proceed to next combination
        BestDrink.erase(currdrink);
    }
}

bool QuxRgb(list<int> &RgbIn)
{
    if(RgbIn.size() <= 1)
        return true;
    
    auto curr = RgbIn.begin();
    auto prev = curr++;
    while(curr != RgbIn.end())
    {
        if(*curr != *prev)
        {
            //if prev and curr are different, replace them with the 3rd val and recursively call again
            int temp_prev_val = *prev;
            int temp_curr_val = *curr;
            int newval = (*curr + *prev) % 4;
            if(!newval)
                newval = 2;
            
            RgbIn.erase(prev);
            *curr = newval;
            
            if( QuxRgb(RgbIn) )
                return true;
            
            //add the deleted values back and proceed to next index
            *curr = temp_curr_val;
            RgbIn.insert(curr, temp_prev_val);
        }
        prev = curr++;
    }
               
    return false;
}

bool FixedPntArray(const vector<int> &In, const int begin, const int end)
{
    if (begin > end)
        return false;
    
    int mid = end + begin;
    mid >>= 1;
    
    if(In[begin] > static_cast<int>(In.size()) || In[end] < 0)
        return false;
    
    if(In[mid] == mid)
        return true;
    if(In[mid] < mid)
        return FixedPntArray(In, mid+1, end);
    return FixedPntArray(In, begin, mid-1);
}

bool MyStringSortComp(const char &A, const char &B)
{
    const string S2 = "cae";
    unordered_map<char, int> sortmap;
    for(int i = 0; i < S2.size(); ++i)
        sortmap.insert(make_pair(S2[i], i));
    
    int Aidx = INT_MAX, Bidx = INT_MAX;
    if(sortmap.find(A) != sortmap.end())
        Aidx = sortmap[A];
    if(sortmap.find(B) != sortmap.end())
        Bidx = sortmap[B];
    return(Aidx < Bidx);
}

int EarliestReachableStep(const vector<int> &MaxReach, const int low, const int high, const int find_idx)
{
    for(int i = low; i <= high; ++i)
        if(MaxReach[i] >= find_idx)
            return i;
    return INT_MAX;
}

int FindQuotient(const int dividend, const int divisor)
{
    if(dividend < divisor)
        return 0;
    
    int powerof2 = 0;
    int divisor_pow = divisor << powerof2;
    while(dividend >= divisor_pow)
    {
        ++powerof2;
        divisor_pow <<= 1;
    }
    powerof2 = 1 << (powerof2-1);
    divisor_pow >>= 1;
    return (powerof2)+FindQuotient(dividend-divisor_pow, divisor);
}

void FindDecodeNumWays(const int Num, int &NumWays)
{
    if(0 == Num)
        return;
    
    int curr_plus_next_digit = Num % 100;
    if(curr_plus_next_digit > 0 && curr_plus_next_digit <= 9)
    {
        ++NumWays;
        FindDecodeNumWays(Num/10, NumWays);
    }
    else if(curr_plus_next_digit > 9 && curr_plus_next_digit <= 26)
    {
        NumWays += 2;
        FindDecodeNumWays(Num/100, NumWays);
    }
    else
        FindDecodeNumWays(Num/10, NumWays);
}

#if 0
    //2 options - include just last digit or last digit + next digit
    int curr_digit, ret1 = INT_MIN;
    curr_digit = Num % 10;
    ret1 = FindDecodeNumWays(Num/10);
    
    int curr_plus_next_digit = INT_MAX, ret2 = INT_MIN;
    if(Num >= 10)
    {
        curr_plus_next_digit = Num % 100;
        if(curr_plus_next_digit <= 26)
            ret2 = FindDecodeNumWays(Num/100);
    }
    
    if(curr_plus_next_digit > 9 && curr_plus_next_digit <= 26 )
        return (max(ret1,ret2) + 2);
    else if (curr_plus_next_digit > 26)
        return (max(ret1,ret2));
    else
        return (max(ret1,ret2) + 1);
}
#endif

bool RearrangeRepeatChar(list<char> &CharList, string &OpStr)
{
    if(CharList.empty())
        return true;
    
    auto size = CharList.size();
    
    //repeat (find different options) for all elements of the list
    for(auto i = 0; i < size; ++i)
    {
        //store the first element of the list and remove it from the list.
        //add it back to the end of the list at the end of recursion
        auto front_item = CharList.front();
        CharList.pop_front();
        
        if(OpStr.empty() || OpStr.back() != front_item)
        {
            OpStr.push_back(front_item);
            if(RearrangeRepeatChar(CharList, OpStr))
                return true;
        }
        
        //add the popped item at the back of the list
        CharList.push_back(front_item);
    }
    
    //pop last output char as this combination does not work
    OpStr.pop_back();
    return false;
}

void FindAllKeypadComb(const int idx, const string &InNum, const vector<string> &NumAlphTable, vector<string> &KeypadOp, string &TempStr)
{
    if(idx >= InNum.size())
    {
        KeypadOp.emplace_back(TempStr);
        return;
    }
    
    //If the digit is 0 or 1, then proceed to the next index
    if('0' == InNum[idx] || '1' == InNum[idx])
        FindAllKeypadComb(idx+1, InNum, NumAlphTable, KeypadOp, TempStr);
    
    for(int i = 0; i < NumAlphTable[InNum[idx]-'0'].size(); ++i)
    {
        TempStr.push_back(NumAlphTable[InNum[idx]-'0'][i]);
        FindAllKeypadComb(idx+1, InNum, NumAlphTable, KeypadOp, TempStr);
        TempStr.pop_back();
    }
}

int CntNumAttempsEggs(const int numEggs, const int numFloors, unordered_map<string,int> &MemoEgg)
{
    if(numEggs == 1)
        return numFloors;
    if(numFloors <= 1)
        return numFloors;
    
    string key = to_string(numEggs) + ":" + to_string(numFloors);
    if(MemoEgg.find(key) != MemoEgg.end())
        return MemoEgg[key];
 
    int ret = INT_MAX, temp;
    for(int i = 1; i <= numFloors; ++i)
    {
        temp = max(CntNumAttempsEggs(numEggs-1, i-1, MemoEgg),
                   CntNumAttempsEggs(numEggs, numFloors-i, MemoEgg));
        ret = min(ret, temp);
    }
    
    MemoEgg.insert(make_pair(key, ret+1));
    
    return ret+1;
}

int GetNumberOfDigits(int A)
{
    int num = 0;
    while(A)
    {
        A /= 10;
        ++num;
    }
    return num;
}

bool MyCompSortForDcp228(int A, int B)
{
    int Adigit = 0, Bdigit = 0;
    int Anumdigits = GetNumberOfDigits(A), Bnumdigits = GetNumberOfDigits(B);
    while(Anumdigits && Bnumdigits)
    {
        --Anumdigits;
        --Bnumdigits;
        Adigit = A / pow(10,Anumdigits);
        Bdigit = B / pow(10,Bnumdigits);
        A -= (Adigit * pow(10,Anumdigits));
        B -= (Bdigit * pow(10,Bnumdigits));
        if(Adigit != Bdigit)
            return (Bdigit < Adigit);
    }
    
    if(Anumdigits == Bnumdigits)
        return true; //this doesnt matter as at this point, both A and B will be the same
    
    if(Anumdigits > Bnumdigits)
    {
        //A has more digits and larger than B. Compare the next digit of A with previous digit of B
        Adigit = A / pow(10,Anumdigits-1);
    }
    else
    {
        //B has more digits and larger than A. Compare the next digit of B with previous digit of A
        Bdigit = B / pow(10,Bnumdigits-1);
    }
    return (Adigit >= Bdigit ? true : false);
}

bool CrackSafeDfs(const int &K, const int &N, unordered_set<string> &Visited, string &OpPwd)
{
    if(Visited.size() >= pow(K,N))
        return true;
    string PwdSubstr = OpPwd.substr(OpPwd.size()-N+1,N-1);
    for(int i = 0; i < K; ++i)
    {
        string NextPwd = PwdSubstr + static_cast<char>('0'+i);
        if(Visited.find(NextPwd) == Visited.end())
        {
            Visited.insert(NextPwd);
            OpPwd.push_back(static_cast<char>('0'+i));
            if( CrackSafeDfs(K,N,Visited,OpPwd) )
                return true;
            else
            {
                Visited.erase(NextPwd);
                OpPwd.pop_back();
            }
        }
    }
    return false;
}

void RunDfsArray(const int idx, vector<bool> &Visited, stack<char> &TopoSortOp, const vector<char> &CharArray)
{
    if(idx >= Visited.size())
        return;
    Visited[idx] = true;
    
    //Goto next even index
    if(idx % 2 == 0)
        if(idx+2 < Visited.size() && Visited[idx+2] == false)
            RunDfsArray(idx+2, Visited, TopoSortOp, CharArray);
    
    //Goto next multiple of 3 index
    if(idx % 3 == 0)
        if(idx+3 < Visited.size() && Visited[idx+3] == false)
            RunDfsArray(idx+3, Visited, TopoSortOp, CharArray);
    
    //Since idx has no more next ngbrs, add to topological stack output
    TopoSortOp.push(CharArray[idx]);
}

string SplitSentance(const string &inString, const unordered_set<string> &WordDict, int idx)
{
    string result;
    
    if(idx >= inString.size())
        return "";
    string temp;
    for(int i = idx; i < inString.size(); ++i)
    {
        //cout << inString.substr(idx,i-idx) << endl;
        temp.push_back(inString[i]);
        if(WordDict.find(temp) != WordDict.end())
        {
            result = SplitSentance(inString, WordDict, i+1);
            if(result != "")
            {
                temp.push_back(' ');
                temp.append(result);
                cout << temp << endl;
                return temp;
            }
        }
    }
    if(WordDict.find(temp) != WordDict.end())
        return temp;
    return "";
}

int FindKthElemToDelete(const int N, const int K)
{
    if(1 == N)
        return 0;
    return( (FindKthElemToDelete(N-1,K) + K)  % N );
}

int MultByPowerOf10(int num, int pow)
{
    int ret = num, temp;
    while(pow)
    {
        temp = (ret << 3) + (ret << 1); //multiply ret by 10
        ret = temp;
        --pow;
    }
    return ret;
}
    
void ComputeMatrixPower(const vector<vector<int>> &A, const int N, vector<vector<int>> &Output)
{
    if(N == 0)
        return;
    if(N == 1)
    {
        Output = A;
        return;
    }
    
    //Compute power N/2 and multiply it twice
    ComputeMatrixPower(A,N/2,Output);
    vector<vector<int>> Temp;
    MultiplyMatrix(Output, Output, Temp);
    
    //if N is odd, multiply input again
    if(N % 2)
        MultiplyMatrix(Temp, A, Output);
    else
        Output = Temp;
}

void MultiplyMatrix(const vector<vector<int>> &A, const vector<vector<int>> &B, vector<vector<int>> &Output)
{
    //Sanity checks - size of all rows of both input matrices should be the same
    //Number of columns of 1st matrix should be same as number of rows of 2nd matrix
    if(A.empty() || B.empty())
        return;
    
    int Arow = A.size();
    int Acol = A[0].size();
    for(const auto &i : A)
        if(i.size() != Acol)
            return;
    
    int Brow = B.size();
    int Bcol = B[0].size();
    for(const auto &i : B)
        if(i.size() != Bcol)
            return;
    
    if(Acol != Brow)
        return;

    Output.resize(Arow);
    for(int i = 0; i < Output.size(); ++i)
        Output[i].resize(Bcol);
    
    for(int i = 0; i < Arow; ++i)
        for(int j = 0; j < Bcol; ++j)
        {
            Output[i][j] = 0;
            for(int k = 0; k < Acol; ++k)
                Output[i][j] += (A[i][k] * B[k][j]);
        }
}

int FindContOnesLocFromLsb(int A)
{
    int ret = INT_MIN;
    int bin_idx = 0, curr_digit, prev_digit;
    prev_digit = A & 1;
    A >>= 1;
    while(A)
    {
        ++bin_idx;
        curr_digit = A & 1;
        A >>= 1;
        
        if(curr_digit & prev_digit)
            //2 back-to-back digits are 1
            ret = bin_idx;
        else
        {
            //2 back-to-back digits are not 1. If the 1st of back-to-back digits with 1 have been identified,
            //return
            if(INT_MIN != ret)
                break;
        }
        prev_digit = curr_digit;
    }
    return ret;
}

shared_ptr<BstType> IsNodePresentBst(const shared_ptr<BstType> Curr, const int A)
{
    if(nullptr == Curr)
        return Curr;
    
    if(A == Curr->key)
        return Curr;
    else if(A < Curr->key)
        return IsNodePresentBst(Curr->left, A);
    else
        return IsNodePresentBst(Curr->right, A);
}

int BstLevelCnt(const shared_ptr<BstType> Curr, const int A)
{
    if(nullptr == Curr)
        return INT_MIN;
    
    if(A == Curr->key)
        return 0;
    else if(A < Curr->key)
        return BstLevelCnt(Curr->left, A) + 1;
    else
        return BstLevelCnt(Curr->right, A) + 1;
}

void FindIpAddr(const string input, const int stringidx, int IpByteVal, vector<int> &IpAddr, int IpBytePos)
{
    static int loopctr = 0;
    
    if(stringidx < 0)
    {
        if(4 == IpAddr.size())
        {
            //valid IP address
            cout << "!!!!!!!!!!!!valid IP!!!!!!!!!!!" << endl;
            PrintDebugIpAddrData(stringidx, IpByteVal, IpAddr, IpBytePos, loopctr);
        }
        else
        {
            cout << "stringidx < 0" << endl;
            PrintDebugIpAddrData(stringidx, IpByteVal, IpAddr, IpBytePos, loopctr);
        }
        return;
    }
    
    if(IpAddr.size() > 4)
    {
        cout << "IpAddr.size() > 4" << endl;
        PrintDebugIpAddrData(stringidx, IpByteVal, IpAddr, IpBytePos, loopctr);
        return;
    }
    
    if(IpAddr.size() == 4 && stringidx >= 2)
    {
        //we have reached size of 4 bytes and still have 3+ digits to process. So its an invalid combination
        cout << "IpAddr.size() == 4 && stringidx >= 2" << endl;
        PrintDebugIpAddrData(stringidx, IpByteVal, IpAddr, IpBytePos, loopctr);
        return;
    }
    
    ++loopctr;
    int digit = static_cast<int>(input[stringidx]) - static_cast<int>('0');
    
    cout << "No errors. Digit " << digit << endl;
    PrintDebugIpAddrData(stringidx, IpByteVal, IpAddr, IpBytePos, loopctr);
    
    //put current digit in IpAddr and advance to next input only if BytePos == 0
    if(0 == IpBytePos)
    {
        IpAddr.push_back(digit);
        FindIpAddr(input, stringidx-1, 0, IpAddr, 0);
        IpAddr.pop_back(); //remove the added digit for next combination
    }
    
    //add current digit to previous byte value and add to IpAddr and advance to next input
    //do this only if current digit is non-zero, and bytepos != 0
    IpByteVal += (pow(10,IpBytePos)*digit);
    if(IpByteVal <= 255 && 0 != digit && 0 != IpBytePos)
    {
        IpAddr.push_back(IpByteVal);
        FindIpAddr(input, stringidx-1, 0, IpAddr, 0);
        IpAddr.pop_back(); //remove the added digit for next combination
    }
    
    //add current digit to previous byte value and continue to next digit withtout adding to IpAddr
    //do this only if we have not reached last digit yet
    if(IpByteVal <= 255 && stringidx > 0)
    {
        FindIpAddr(input, stringidx-1, IpByteVal, IpAddr, IpBytePos+1);
    }
}

void PrintDebugIpAddrData(const int &stringidx, const int &IpByteVal, const vector<int> &IpAddr, const int &IpBytePos, const int &LoopCtr)
{
    cout << "Idx " << stringidx << " IpByteVal " << IpByteVal << " BytePos " << IpBytePos << " loopctr " << LoopCtr << endl;
    for(const auto &i:IpAddr)
        cout << i << ":";
    cout << endl << endl;
}

int FindCollatzCnt(const int N, unordered_map<int,int> &CollatzMap)
{
    if(CollatzMap.find(N) != CollatzMap.end())
        return CollatzMap[N];
        
    int RetCnt = 1, temp = N;
    while(temp > 1)
    {
        if(temp%2 == 0)
            temp = temp/2;
        else
            temp = (3*temp)+1;
        
        if(CollatzMap.find(temp) != CollatzMap.end())
        {
            RetCnt += CollatzMap[temp];
            break;
        }
        else
            ++RetCnt;
    }
    
    CollatzMap.insert(make_pair(N, RetCnt));
    return RetCnt;
}

void PrintLinkedList(const shared_ptr<LinkListNodeType> Head)
{
    shared_ptr<LinkListNodeType> Temp = Head;
    cout << "Print LL ";
    while(Temp)
    {
        cout << Temp->val << " ";
        Temp = Temp->next;
    }
    cout << endl;
}

shared_ptr<LinkListNodeType> AddNodeToLlHead(shared_ptr<LinkListNodeType> Head, shared_ptr<LinkListNodeType> Curr)
{
    if(nullptr == Curr)
        return Head;
    Curr->next = Head;
    return Curr;
}

shared_ptr<LinkListNodeType> AddNodeToLlTail(shared_ptr<LinkListNodeType> Tail, shared_ptr<LinkListNodeType> Curr)
{
    if(nullptr == Curr)
        return Tail;
    Curr->next = nullptr;
    if(Tail)
        Tail->next = Curr;
    return Curr;
}

int GetDigit(const int &Input, const int idx)
{
    int digit = Input/static_cast<int>(pow(10,idx));
    digit %= 10;
    return digit;
}

int FindMinElem(const vector<int> &Input, const int start, const int end)
{
    cout << start << " " << end << endl;
    
    if(start > end)
        return INT_MAX;
    if(start == end)
        return Input[start];
    
    int middle = (start+end)/2;
    
    if(Input[middle] >= Input[start] && Input[end] >= Input[middle+1])
        return min(Input[start], Input[middle+1]);
    else if(Input[middle] < Input[start])
        return FindMinElem(Input, start, middle);
    else
        return FindMinElem(Input, middle+1, end);
}

bool TriangleErrorCheck(const vector<vector<int>> &Input)
{
    for(int i = 0; i < Input.size(); ++i)
        if(Input[i].size() != i+1)
            return false;
    return true;
}

pair<int,int> TriangleGetRight(const pair<int,int> Curr, const int M)
{
    if(Curr.first >= M-1)
        return make_pair(INT_MAX, INT_MAX);
    return make_pair(Curr.first+1, Curr.second+1);
}

pair<int,int> TriangleGetLeft(const pair<int,int> Curr, const int M)
{
    if(Curr.first >= M-1)
        return make_pair(INT_MAX, INT_MAX);
    return make_pair(Curr.first+1, Curr.second);
}

int TriangleMaxSum(const vector<vector<int>> &Input, const pair<int,int> Curr)
{
    if(Curr.first == INT_MAX || Curr.second == INT_MAX)
        return 0;
    pair<int,int> Left = TriangleGetLeft(Curr, Input.size());
    pair<int,int> Right = TriangleGetRight(Curr, Input.size());
    
    return (Input[Curr.first][Curr.second] + max( TriangleMaxSum(Input,Left), TriangleMaxSum(Input,Right) ));
}

void ReverseArray(vector<int> &Input, const int start, const int end)
{
    if(start < 0 || start >= Input.size() || end < 0 || end >= Input.size())
        return;
    
    int temp, half_size = (end-start)/2;
    for(int i = 0; i <= half_size; ++i)
    {
        temp = Input[start+i];
        Input[start+i] = Input[end-i];
        Input[end-i] = temp;
    }
}

int SumOfNodesRepetitions(const shared_ptr<BstType> Curr, unordered_map<int, int> &SumCntMap)
{
    if(nullptr == Curr)
        return 0;
    int Sum = SumOfNodesRepetitions(Curr->left, SumCntMap) + SumOfNodesRepetitions(Curr->right, SumCntMap) + Curr->key;
    SumCntMap[Sum]++;
    return(Sum);
}

typedef pair<int,pair<int, int>> KTH_SMALL_TYPE;
struct CompareKthSmall
{
    bool operator()(const KTH_SMALL_TYPE &a, const KTH_SMALL_TYPE &b)
    {
        return a.first > b.first;
    }
};

int KthSmallest2d(const vector<vector<int>> &arr, const int k)
{
    int M = arr.size(), N = arr[0].size();
    
    //priority_queue<KTH_SMALL_TYPE, vector<KTH_SMALL_TYPE>, CompareKthSmall> pq;
    priority_queue<KTH_SMALL_TYPE, vector<KTH_SMALL_TYPE>, greater<KTH_SMALL_TYPE>> pq;
    
    //add 1st column in pq
    for(int i = 0; i < M; i++)
        pq.push( make_pair(arr[i][0], make_pair(i,0)) );
    
    int x = k, ans = 0;
    while(x--)
    {
        int e = pq.top().first;
        int i = pq.top().second.first;
        int j = pq.top().second.second;
        ans = e;
        pq.pop();
        if(j < N-1)
            pq.push( make_pair( arr[i][j+1], make_pair(i,j+1) ) );
    }
    return ans;
}

int FindIdxSortedArray(const vector<vector<int>> &Input, const int FindNum, const int row, const int left, const int right)
{
    if(FindNum <= Input[row][left])
        return left;
    if(FindNum > Input[row][right])
        return right+1;
    if(left >= right)
        return right;
    
    int middle = (left + right)/2;
    
    if(FindNum == Input[row][middle])
       return middle;
    else if(FindNum > Input[row][middle])
        return FindIdxSortedArray(Input, FindNum, row, middle+1, right);
    else
        return FindIdxSortedArray(Input, FindNum, row, left, middle-1);
}

int ComputeInversions(const vector<pair<int,int>> &LineSegments)
{
    vector<int> TempMerge (LineSegments.size(), 0), MergetInput;
    for(const auto &i : LineSegments)
        MergetInput.push_back(i.second);
    return MergeSortInvCnt(MergetInput, TempMerge, 0, LineSegments.size()-1);
}

int MergeSortInvCnt(vector<int> &Input, vector<int> &TempMerge, const int LeftStart, const int RightEnd)
{
    if(LeftStart >= RightEnd)
        return 0;
    
    int Ret = 0, Middle = (RightEnd+LeftStart)/2;
    Ret += MergeSortInvCnt(Input, TempMerge, LeftStart, Middle);
    Ret += MergeSortInvCnt(Input, TempMerge, Middle+1, RightEnd);
    Ret += MergeSortInvCntCombine(Input, TempMerge, LeftStart, RightEnd);
    return Ret;
}

int MergeSortInvCntCombine(vector<int> &Input, vector<int> &TempMerge, const int LeftStart, const int RightEnd)
{
    if(LeftStart >= RightEnd)
        return 0;
    
    int LeftEnd = (RightEnd+LeftStart)/2;
    int RightStart = LeftEnd+1;
    int leftidx = LeftStart, rightidx = RightStart, combine_idx = LeftStart;
    int TotalInv = 0;
    
    //combine left and right arrays
    while(leftidx <= LeftEnd && rightidx <= RightEnd)
    {
        if(Input[leftidx] <= Input[rightidx])
        {
            TempMerge[combine_idx] = Input[leftidx];
            ++leftidx;
        }
        else
        {
            TempMerge[combine_idx] = Input[rightidx];
            ++rightidx;
            TotalInv += (LeftEnd - leftidx + 1);
        }
        ++combine_idx;
    }
    
    //append any leftover rightarray elements
    while(rightidx <= RightEnd)
    {
        TempMerge[combine_idx] = Input[rightidx];
        ++rightidx;
        ++combine_idx;
    }
    
    //append any leftover leftarray elements and add inversions
    while(leftidx <= LeftEnd)
    {
        TempMerge[combine_idx] = Input[leftidx];
        ++leftidx;
        ++combine_idx;
        //TotalInv += (RightEnd - RightStart + 1);
    }
    
    //copy merged output
    for(int i = LeftStart; i <= RightEnd; ++i)
        Input[i] = TempMerge[i];
    
    return TotalInv;
}

bool IsArrayEnd(const vector<int> &Input, const int curridx, unordered_map<int, bool> &Memo)
{
    if(curridx < 0 || curridx >= Input.size())
        return false;
    if(curridx == Input.size()-1)
        return true;
    if(0 == Input[curridx])
        return false;
    
    bool Ret = false;
    for(int i = 1; i <= Input[curridx]; ++i)
    {
        if(Memo.find(curridx) != Memo.end())
            return Memo[curridx];
        if(IsArrayEnd(Input, curridx+i, Memo))
        {
            Ret = true;
            break;
        }
    }
    Memo.insert(make_pair(curridx, Ret));
    return Ret;
}

int MaxSumSubarray(const vector<int> &Input, const int start, const int end)
{
    if(start < 0 || end < 0 || start >= Input.size() || end >= Input.size())
        return INT_MIN;
    if(start == end)
        return Input[start];
    
    int MaxHere, MaxSoFar;
    MaxHere = MaxSoFar = Input[start];
    int i = start + 1;
    
    while(1)
    {
        MaxHere = max(MaxHere + Input[i], Input[i]);
        MaxSoFar = max(MaxSoFar, MaxHere);
        if(i == end)
            break;
        ++i;
        i %= Input.size();
    }
    
    return MaxSoFar;
}

int NumSqInts(const int SqSumNum, const vector<int> &NumArray, int NumInts, int MaxIdx)
{
    if(SqSumNum == 0)
        return NumInts;
    if(SqSumNum < 0 || MaxIdx < 0)
        return INT_MAX;
    
    //Memoization
    static unordered_map<string,int> Memo;
    string MemoKey = to_string(SqSumNum) + ":" + to_string(MaxIdx);
    if(Memo.find(MemoKey) != Memo.end())
        return Memo[MemoKey];

    static int LoopCtr = 0;
    cout << "LoopCtr " << ++LoopCtr << " Sum " << SqSumNum << " MaxIdx " << MaxIdx << endl;
    
    int Ret = INT_MAX, temp;
    for(int i = MaxIdx; i >= 0 ; --i)
    {
        if(SqSumNum < NumArray[i]*NumArray[i])
            temp = NumSqInts(SqSumNum, NumArray, NumInts, i-1);
        else
            temp = min( NumSqInts(SqSumNum-(NumArray[i]*NumArray[i]), NumArray, NumInts+1, MaxIdx),
                       NumSqInts(SqSumNum, NumArray, NumInts, MaxIdx-1) );
        Ret = min(Ret, temp);
    }
    
    Memo.insert(make_pair(MemoKey, Ret));
    
    return Ret;
}

int LenLongestDistint(const vector<int> &Input)
{
    int start_idx = 0, end_idx = 0;
    int MaxDistLen = 0;
    unordered_map<int,int> DistinctSet;
    while(start_idx <= end_idx && start_idx < Input.size() && end_idx < Input.size())
    {
        auto Itr = DistinctSet.find(Input[end_idx]);
        if(Itr == DistinctSet.end())
        {
            DistinctSet.insert(make_pair(Input[end_idx],end_idx));
            ++end_idx;
        }
        else
        {
            MaxDistLen = max( MaxDistLen, static_cast<int>(DistinctSet.size()) );
            
            //Remove all elements upto and including the repeating element, and start the start_idx from the next idx after
            // the 1st occurence of the repeating element
            start_idx = Itr->second + 1;
            vector<unordered_map<int,int>::iterator> EraseMapItrVec;
            for(auto i = DistinctSet.begin(); i != DistinctSet.end(); ++i)
            {
                cout << i->second << " " << i->first << endl;
                if(i->second < start_idx)
                {
                    cout << "erase " << i->second << " " << i->first << endl;
                    EraseMapItrVec.push_back(i);
                }
            }
            for(auto &i : EraseMapItrVec)
                DistinctSet.erase(i);
        }
    }
    
    return MaxDistLen;
}

int MinPalinPartCnt(const string &Input, unordered_map<string,int> &Memo, const int i, const int j)
{
    if(i == j)
        return 0;
    if(i > j || i >= Input.size() || j >= Input.size())
        return INT_MAX;
    
    //Memoization
    string MemoKey = to_string(i) + ":" + to_string(j);
    if(Memo.find(MemoKey) != Memo.end())
        return Memo[MemoKey];
    
    int MinPart = INT_MAX, FirstPart = 0, SecondPart = 0;
    static int LoopCnt = 0;
    cout << "LoopCnt " << ++LoopCnt << " i " << i << " j " << j << endl;
    
    if(IsPalindrome(Input,i,j))
        MinPart = 0;
    else
    {
        for(int k = i; k < j; ++k)
        {
            int FirstPart = MinPalinPartCnt(Input, Memo, i, k);
            int SecondPart = MinPalinPartCnt(Input, Memo, k+1, j);
            if(INT_MAX == FirstPart || INT_MAX == SecondPart)
                continue;
            else
                MinPart = min(MinPart, FirstPart+SecondPart+1);
        }
    }
    
    Memo.insert(make_pair(MemoKey, MinPart));
    return MinPart;
}

bool IsPalindrome(const string &Input, const int i, const int j)
{
    if(i > j || i >= Input.size() || j >= Input.size())
        return false;
    if(i == j)
        return true;
    for(int k = 0; k <= (j-i)/2; ++k)
        if(Input[k+i] != Input[j-k])
            return false;
    return true;
}

int MinSubsetFn(const int TotalSum, const int SubarraySum, int idx, const vector<int> &Input)
{
    static int LoopCnt = 0;
    
    if(idx < 0)
        return abs(TotalSum-SubarraySum-SubarraySum);
    
    cout << " SubarraySum " << SubarraySum << " idx " << idx << " LoopCtr " << ++LoopCnt << endl;
    
    return(min( MinSubsetFn(TotalSum, SubarraySum+Input[idx], idx-1, Input),
               MinSubsetFn(TotalSum, SubarraySum, idx-1, Input) ));
}

void GenerateGrayCode(const int num, const int Len, unordered_set<int> &GrayCodeSet)
{
    static int Loopcnt = 0;
    cout << "GenerateGrayCode: Loopcnt " << ++Loopcnt << " Num: " << num << endl;
    
    int mask, temp;
    for(int i = 0; i < Len; ++i)
    {
        mask = 1 << i; //mask to flip 1 bit
        temp = num ^ mask; //flip a bit
        if(GrayCodeSet.find(temp) == GrayCodeSet.end())
        {
            //temp not present in set. Add to set and flip starting from 0th bit
            GrayCodeSet.insert(temp);
            GenerateGrayCode(temp, Len, GrayCodeSet);
        }
    }
}

bool IsSubtreeCheck(const shared_ptr<BstType> CurrT, const shared_ptr<BstType> CurrS, const shared_ptr<BstType> &BstHeadS )
{
    bool a = false, b = false;
    
    if(nullptr == CurrT)
        return (nullptr == CurrS ? true : false);
    if(nullptr == CurrS)
        return false;
    
    if(CurrS->key != CurrT->key)
    {
        a = IsSubtreeCheck(CurrT->left, BstHeadS, BstHeadS);
        if(!a)
            b = IsSubtreeCheck(CurrT->right, BstHeadS, BstHeadS);
        return (a || b);
    }
    else
    {
        a = IsSubtreeCheck(CurrT->left, CurrS->left, BstHeadS);
        b = IsSubtreeCheck(CurrT->right, CurrS->right, BstHeadS);
        return (a && b);
    }
}

pair<int,int> TreeLevelMinSum(const shared_ptr<BstType> Curr, int sum, int level)
{
    if(nullptr == Curr)
        return(make_pair(INT_MAX, INT_MAX));
    
    sum += Curr->key;
    ++level;
    
    if(nullptr == Curr->left && nullptr == Curr->right)
        return(make_pair(sum, level));
    
    auto right = TreeLevelMinSum(Curr->right, sum, level);
    auto left = TreeLevelMinSum(Curr->left, sum, level);
    
    return(right.first < left.first ? right : left);
}

void MergeSortedVecs(const vector<int> &A, const vector<int> &B, vector<int> &Op)
{
    int Aidx = 0, Bidx = 0, Opidx = 0;
    if(0 == A.size())
    {
        Op = B;
        return;
    }
    if(0 == B.size())
    {
        Op = A;
        return;
    }
    Op.resize(A.size()+B.size(), 0);
    while(Opidx < Op.size())
    {
        Op[Opidx++] = (A[Aidx] <= B[Bidx] ? A[Aidx++] : B[Bidx++]);
        if(Aidx >= A.size() || Bidx >= B.size())
            break;
    }
    while(Aidx < A.size())
        Op[Opidx++] = A[Aidx++];
    while(Bidx < B.size())
        Op[Opidx++] = B[Bidx++];
}

void ReverseWord(string &word, int begin, int end)
{
    char temp;
    while(begin < end)
    {
        temp = word[begin];
        word[begin] = word[end];
        word[end] = temp;
        ++begin;
        --end;
    }
}

bool HopFunction(const int curridx, const vector<int> &HopInput, unordered_map<int,bool> &HopMemo)
{
    static int loopcnt = 0;
    ++loopcnt;
    cout << curridx << " " << loopcnt << endl;
    
    bool ret = false;
    if(curridx >= HopInput.size())
        return false;
    if(HopInput.size()-1 == curridx)
        return true;
    if(0 == HopInput[curridx])
        return false;
    if(HopMemo.find(curridx) != HopMemo.end())
        return HopMemo[curridx];
    for(int i = 1; i <= HopInput[curridx]; ++i)
    {
        ret = HopFunction(curridx+i, HopInput, HopMemo);
        if(ret)
            break;
    }
    HopMemo.insert(make_pair(curridx, ret));
    return ret;
}

void AllPermOfNumVec(vector<vector<int>> &AllPermVecOp, vector<int> DigitCnt, int level, vector<int> TempVecOp)
{
    if(level < 0)
    {
        AllPermVecOp.push_back(TempVecOp);
        return;
    }
    for(int i = 0; i < 128; ++i)
        if(DigitCnt[i]) //if there is a non-zero cnt value in digit cnt array, add to output vector
        {
            TempVecOp.push_back(i-'0');
            --DigitCnt[i]; //add to output vector and decrement the cnt, and proceed to next digit
            AllPermOfNumVec(AllPermVecOp,DigitCnt,level-1,TempVecOp);
            ++DigitCnt[i];
            TempVecOp.pop_back(); //increment the cnt and remove the val from output variable for next combination
        }
}

//Binary function that accepts two elements in the range as arguments, and returns a value convertible to bool.
//The value returned indicates whether the element passed as first argument is considered to go before the second
//in the specific strict weak ordering it defines.
bool MyCompDecreasingOrder(const int i, const int j)
{
    return (j<i);
}

shared_ptr<BstType> LcaBst(const shared_ptr<BstType> Curr, const int A, const int B)
{
    if(nullptr == Curr)
        return nullptr;
    if(max(A,B) < Curr->key)
        return LcaBst(Curr->left, A, B);
    else if(min(A,B) > Curr->key)
        return LcaBst(Curr->right, A, B);
    else
        return Curr;
}

shared_ptr<BstType> LcaBt(const shared_ptr<BstType> Curr, const int A, const int B)
{
    if(nullptr == Curr)
        return nullptr;
    if(Curr->key == A || Curr->key == B)
        return Curr;
    shared_ptr<BstType> left = LcaBt(Curr->left, A, B);
    shared_ptr<BstType> right = LcaBt(Curr->right, A, B);
    if(left==nullptr && right==nullptr)
        return nullptr;
    if(left==nullptr)
        return right;
    if(right==nullptr)
        return left;
    return Curr;
}

pair<shared_ptr<BstType>,int> SumLeafToRoot(const shared_ptr<BstType> Curr, const int key)
{
    if(nullptr == Curr)
        return make_pair(nullptr,0);
    if(key == Curr->key)
        return make_pair(Curr,Curr->key);
    
    pair<shared_ptr<BstType>,int> left = SumLeafToRoot(Curr->left, key);
    pair<shared_ptr<BstType>,int> right = SumLeafToRoot(Curr->right, key);
    
    if(left.first != nullptr)
        return make_pair(left.first, left.second+Curr->key);
    else if(right.first != nullptr)
        return make_pair(right.first, right.second+Curr->key);
    else
        return make_pair(nullptr,0);
}

int SumOfAllNodes(const shared_ptr<BstType> Curr)
{
    if(nullptr == Curr)
        return 0;
    return(SumOfAllNodes(Curr->left) + SumOfAllNodes(Curr->right) + Curr->key);
}

int MaxSumAnyLeafToRoot(const shared_ptr<BstType> Curr)
{
    if(nullptr == Curr)
        return 0;
    return( max(MaxSumAnyLeafToRoot(Curr->left), MaxSumAnyLeafToRoot(Curr->right)) + Curr->key );
}

int HeightBinTree(const shared_ptr<BstType> Curr)
{
    if(nullptr == Curr)
        return 0;
    return( max(HeightBinTree(Curr->left), HeightBinTree(Curr->right)) + 1);
}

//O(N) runtime, O(h) space - better for balanced tree (h is height - stored in recursive stack)
void InvertBinaryTreeRecur(shared_ptr<BstType> Curr)
{
    if(nullptr == Curr)
        return;
    
    //Swap left and right tree
    shared_ptr<BstType> Temp = Curr->right;
    Curr->right = Curr->left;
    Curr->left = Temp;
    
    //Do the same for the child nodes
    InvertBinaryTreeRecur(Curr->left);
    InvertBinaryTreeRecur(Curr->right);
}

//O(N) runtime, O(N/2) space - better for unbalanced tree
void InvertBinaryTreeIter(shared_ptr<BstType> Head)
{
    if(nullptr == Head)
        return;
    
    queue<shared_ptr<BstType>> PtrQ;
    PtrQ.push(Head);
    
    while(!PtrQ.empty())
    {
        shared_ptr<BstType> Curr = PtrQ.front(), Temp;
        PtrQ.pop();
        
        //Swap left and right tree
        Temp = Curr->right;
        Curr->right = Curr->left;
        Curr->left = Temp;
        
        //Do the same for the child nodes
        if(Curr->left)
            PtrQ.push(Curr->left);
        if(Curr->right)
            PtrQ.push(Curr->right);
    }
}

pair<shared_ptr<BstType>,int> FindDeepestNodeTree(const shared_ptr<BstType> Curr, const int level)
{
    if(nullptr == Curr)
        return make_pair(nullptr,0);
    if(nullptr == Curr->left && nullptr == Curr->right)
        return make_pair(Curr,level);
    
    pair<shared_ptr<BstType>,int> left = FindDeepestNodeTree(Curr->left, level+1);
    pair<shared_ptr<BstType>,int> right = FindDeepestNodeTree(Curr->right, level+1);
    
    return (left.second > right.second ? left : right );
}

pair<bool,int> FindLargestBst(const shared_ptr<BstType> Curr, const int Min, const int Max)
{
    if(nullptr == Curr)
        return make_pair(true,0);
    
    pair<bool,int> left = FindLargestBst(Curr->left, Min, Curr->key);
    pair<bool,int> right = FindLargestBst(Curr->right, Curr->key, Max);
    
    if(Curr->key < Min || Curr->key > Max || !left.first || !right.first)
        return make_pair( false, max(left.second,right.second) );
    else
        return make_pair( true, left.second + right.second + 1 );
}

#if 0
//Insert key to non-binary search tree
void InsertKeyToAnyTree(const int key, const string value)
{
    if(nullptr == BstHead)
    {
        BstHead = make_shared<BstType>();
        BstHead->key = key;
        BstHead->value = value;
        cout << "Insert non-BST head complete: Key: " << key << " Value: " << value << endl;
        return;
    }
    
    queue<shared_ptr<BstType>> PtrQ;
    PtrQ.push(BstHead);
    
    while(!PtrQ.empty())
    {
        shared_ptr<BstType> Curr = PtrQ.front();
        PtrQ.pop();
        
        if(nullptr == Curr->left)
        {
            Curr->left = make_shared<BstType>();
            Curr->left->key = key;
            Curr->left->value = value;
            Curr->size = GetSize(Curr->left) + GetSize(Curr->right) + 1;
            break;
        }
        else if(nullptr == Curr->right)
        {
            Curr->right = make_shared<BstType>();
            Curr->right->key = key;
            Curr->right->value = value;
            Curr->size = GetSize(Curr->left) + GetSize(Curr->right) + 1;
            break;
        }
        else
        {
            PtrQ.push(Curr->left);
            PtrQ.push(Curr->right);
        }
    }
    
    cout << "Insert non-BST complete: Key: " << key << " Value: " << value << endl;
}
#endif //#if 0

//Insert key to non-binary search tree
void InsertKeyToAnyTree(const int key, const string value)
{
    BstHead = InsertKeyToAnyTree(BstHead, key, value);
    cout << "Insert complete: Key: " << key << " Value: " << value << endl;
}

shared_ptr<BstType> InsertKeyToAnyTree(shared_ptr<BstType> Node, const int key, const string value)
{
    //Key is not present. Create a new node, and update size
    if (nullptr == Node)
    {
        shared_ptr<BstType> NewNode = make_shared<BstType>();
        NewNode->key = key;
        NewNode->value = value;
        return NewNode;
    }
    
    int leftsize = GetSize(Node->left);
    int rightsize = GetSize(Node->right);
    
    //To make it balanced, insert in the subtree that has fewer elements
    if (leftsize > rightsize)
        Node->right = InsertKeyToAnyTree(Node->right, key, value);
    else
        Node->left = InsertKeyToAnyTree(Node->left, key, value);
    
    Node->size++;
    return Node;
}

bool IsBst(const shared_ptr<BstType> Curr, const int Min, const int Max)
{
    if(nullptr == Curr)
        return true;
    if(Curr->key < Min || Curr->key > Max)
        return false;
    return( IsBst(Curr->left, Min, Curr->key) && IsBst(Curr->right, Curr->key, Max) );
}

int GetNonuniformRandIdx(const int R, const vector<int> &CumProb, int low, int hi)
{
    int mid;
    
    while (low <= hi)
    {
        mid = (low+hi) >> 1;
        if((0 == mid) && (R <= CumProb[mid]))
            break;
        if(( R > CumProb[mid-1] ) && (R <= CumProb[mid] ))
            break;
        (R <= CumProb[mid]) ? (hi = mid - 1) : (low = mid + 1);
    }
    
    //cout << "Dbg R " << R << " Idx " << mid << endl;
    return mid;
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
    
bool PairSortComp(const pair<int,int> i, const pair<int,int> j)
{
    if(i.first == j.first)
        return (j.second < i.second);
    return (i.first < j.first);
}

bool PairSortComp2(const pair<int,int> i, const pair<int,int> j)
{
    if(i.first == j.first)
        return (i.second < j.second);
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

void InsertKeyToBst(const int key, const string value)
{
    BstHead = InsertKeyToBst(BstHead, key, value);
    cout << "Insert complete: Key: " << key << " Value: " << value << endl;
}

shared_ptr<BstType> InsertKeyToBst(shared_ptr<BstType> Node, const int key, const string value)
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
        Node->right = InsertKeyToBst(Node->right, key, value);
    else
        Node->left = InsertKeyToBst(Node->left, key, value);
    
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
