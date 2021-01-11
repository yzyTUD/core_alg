#pragma once

#include<iostream>
#include<queue>
#include<stack>
#include<cstring>
#include<string>


#define OPENMESH_SUPPORT 0

/*
    using a tree structure.
*/ 
struct Trunk
{
    Trunk *prev;
    std::string str;
    Trunk(Trunk *prev, std::string str)
    {
        this->prev = prev;
        this->str = str;
    }
};
class TreeNode{  
    int v;
    TreeNode* left;
    TreeNode* right;
public:
    TreeNode(int value){
        v=value;
        left=right=NULL;
    }
    void setLeftNode(TreeNode* l){
        left=l;
    }
    void setRightNode(TreeNode* r){
        right=r;
    }
    TreeNode* getLeftNode(){
        return left;
    }
    TreeNode* getRightNode(){
        return right;
    }
    int getValue(){
        return v;
    }
};
class zrTree{ 
    TreeNode* root;
    public:
        void buildTreeFromString(char* ss); 
        void showTrunks(Trunk *p);
        void printTree(TreeNode* TreeNode,Trunk *prev, bool isLeft);
        void firstOrderVisit(TreeNode* TreeNode);
        void inOrderVisit(TreeNode* TreeNode);
        void gen_test_data_and_test(){
            char* sstr = (char*)"(324,)(22,L) (3,R) (4,LL)(281,LLL)(5,LR)()"; 
                //-> warning-free - cast to char* var.
            buildTreeFromString(sstr);
            printTree(root,NULL,false);
            std::cout<<"\n";
            std::cout<<"First Order Visit: \t";firstOrderVisit(root);
            std::cout<<"\nIn Order Visit: \t";inOrderVisit(root);
            std::cout<<"\n";
        }
};
void zrTree::buildTreeFromString(char* ss){
    //evaluate the expression, between ( and ,
    char* sp = std::strchr(ss,'(')+1;
    char* ep = std::strchr(sp,',');
    while(sp!=NULL&& ep!=NULL){
        int n=ep-sp;
        int value=0;
        for(int i=0;i<n;i++){
            int factor=1;
            for(int j=n-1-i;j>0;j--){
                factor *= 10;
            }
            value += (*(sp+i)-'0')*factor;
        }
        //std::cout<<value<<std::endl;
        //extract the lev information 
        sp = ep+1;//, 
        ep = std::strchr(sp,')');
        n=ep-sp;
        //std::cout<<n;
        TreeNode* curNode;
        TreeNode* curParent;
        curParent=root;

        //navigate to :
        for(int i=0;i<n-1;i++){
            if(*(sp+i)=='L'){
                //printf("left");
                curParent=curParent->getLeftNode();
            }else if(*(sp+i)=='R'){
                //right
                curParent=curParent->getRightNode();
            }
        }
        // add TreeNode
        if(0==n){
            //root 
            root=new TreeNode(value);
        }else{
            curNode = new TreeNode(value);
            if(*(sp+n-1)=='L'){ 
                curParent->setLeftNode(curNode);
            }else{
                curParent->setRightNode(curNode);
            }
        }

        sp = std::strchr(sp,'(')+1; 
        ep = std::strchr(sp,',');
    }
}
// Helper function to print branches of the binary tree
void zrTree::showTrunks(Trunk *p)
{
    if (p == NULL)
        return;

    showTrunks(p->prev);

    std::cout << p->str;
}
// Recursive function to print binary tree
// It uses inorder traversal
void zrTree::printTree(TreeNode* TreeNode,Trunk *prev, bool isLeft)
{
    if (TreeNode == NULL)
        return;
    
    std::string prev_str = "    ";
    Trunk *trunk = new Trunk(prev, prev_str);

    printTree(TreeNode->getRightNode(), trunk, true);

    if (!prev)
        trunk->str = "---";
    else if (isLeft)
    {
        trunk->str = ".---";
        prev_str = "   |";
    }
    else
    {
        trunk->str = "`---";
        prev->str = prev_str;
    }

    showTrunks(trunk);
    std::cout << TreeNode->getValue() << std::endl;

    if (prev)
        prev->str = prev_str;
    trunk->str = "   |";

    printTree(TreeNode->getLeftNode(), trunk, false);
}
//some operations to the tree 
void zrTree::firstOrderVisit(TreeNode* TreeNode){
    if(TreeNode==NULL){
        return;
    }
    std::cout<<TreeNode->getValue()<<" ";
    firstOrderVisit(TreeNode->getLeftNode());
    firstOrderVisit(TreeNode->getRightNode());
}
void zrTree::inOrderVisit(TreeNode* TreeNode){
    if(TreeNode==NULL){
        return;
    }
    inOrderVisit(TreeNode->getLeftNode());
    std::cout<<TreeNode->getValue()<<" ";
    inOrderVisit(TreeNode->getRightNode());
}


/* using Graph
edges: node0 node1 val
1 2 10
1 4 5
2 3 5
4 3 8
2 5 5

-->1,2,10,1,4,5,2,3,5,4,3,8,2,5,5,2,7,6

https://www.jianshu.com/p/11e3db610c4b
0,1,1,
0,2,2,
1,2,4,
2,0,3,
2,3,6,
3,3,2

--> 0,1,1,0,2,2,1,2,4,2,0,3,2,3,6,3,3,2

*/
const int MAX_GRAPH_NODES=100;
class Edge{
    public:
        int val;
        int vexid_linked_to;
        Edge* next_edge;

        Edge(int v,int vexid){
            val = v;
            vexid_linked_to=vexid;
            next_edge=NULL;
        }
};
class zrGraph{
    private:
        Edge* edgeList[MAX_GRAPH_NODES];
    public:
        void addNodeAtTail(int source_node,Edge* edge_tobeadded){
            Edge* node_ptr = edgeList[source_node];// get the header ptr 
            if(node_ptr==NULL){
                edgeList[source_node]=edge_tobeadded;
            }else{
                while(node_ptr->next_edge!=NULL){// move to the end of the list 
                    node_ptr=node_ptr->next_edge;
                }
                node_ptr->next_edge=edge_tobeadded;// add the new node as next 
            }
        }
        void build_graph_from_array(int* array,int size);
        void showListedGraph();
        void dfs(int begin_vexid);
        void bfs(int begin_vexid);
        void gen_test_data_and_test(){
            int arr[]={0,1,1,0,2,2,1,2,4,2,0,3,2,3,6,3,3,2};
            build_graph_from_array(arr,sizeof(arr)/sizeof(arr[0]));
            showListedGraph();
            bfs(2);
            dfs(2);
        }
};
void zrGraph::build_graph_from_array(int* array,int size){
    for(int i=0;i<size;i+=3){
        this->addNodeAtTail(array[i],new Edge(array[i+2],array[i+1]));
    }
}
void zrGraph::showListedGraph(){
    for(int i=0;i<10;i++){// as an example, just show 10 lines 
        Edge* tmpEdge = this->edgeList[i];
        if(tmpEdge!=NULL){
            std::cout<<"("<<i<<")";
            while(tmpEdge){
                std::cout<<"->(f:"<<i<<",t:"<<tmpEdge->vexid_linked_to
                    <<",val:"<<tmpEdge->val<<")";
                tmpEdge=tmpEdge->next_edge;
            }
            std::cout<<"\n";
        }
    }
}
/*
    push the beginning item to queue
    while queue not empty:
        visit item by pop() and set to visited
        push the not visited adjacent items of it 
    
    f::use queue as an interface to accessing the items 
    f::do not forget to use the array: visited to control the workflow
        otherweise, just think about two lists, will loop forever 
    f::do not have recur. version!
*/
void zrGraph::bfs(int begin_vexid){
    bool* visited=new bool[MAX_GRAPH_NODES];
    for(int i=0;i<MAX_GRAPH_NODES;i++)
        visited[i]=false;
    std::queue<int> q;
    q.push(begin_vexid);
    //visited[begin_vexid]=true;
    while(!q.empty()){
        //visiting item i and dequeue
        int item=q.front();
        q.pop();
        visited[item]=true;
        std::cout<<item<<"->";//<<std::endl;

        //push the adjcent items to the queue
        Edge* curEdge = this->edgeList[item];
        while(curEdge!=NULL){
            if(visited[curEdge->vexid_linked_to]==false)
                q.push(curEdge->vexid_linked_to);
            curEdge=curEdge->next_edge;
        }
    }std::cout<<"\n";
}
/*
    non-recur. 
*/
void zrGraph::dfs(int begin_vexid){
    bool* visited=new bool[MAX_GRAPH_NODES];
    for(int i=0;i<MAX_GRAPH_NODES;i++)
        visited[i]=false;
    std::stack<int> s;
    s.push(begin_vexid);
    //visited[begin_vexid]=true;
    while(!s.empty()){
        //visiting item i and dequeue
        int item=s.top();
        s.pop();
        visited[item]=true;
        std::cout<<item<<"->";//<<std::endl;

        /*
            push the adjcent items to the queue - last item first

        */
        // Edge* curEdge = this->edgeList[item];
        // while(curEdge!=NULL){
        //     if(visited[curEdge->vexid_linked_to]==false)
        //         s.push(curEdge->vexid_linked_to);
        //     curEdge=curEdge->next_edge;
        // }

        /*
            push the adjcent items to the queue - last item last
                push in another direction, still not the same as queue
                f::inverse visit of a linked list item use stack
        */
        std::stack<int> helperstack;
        int helperval;
        Edge* curEdge = this->edgeList[item];
        while(curEdge!=NULL){
            if(visited[curEdge->vexid_linked_to]==false)
                helperstack.push(curEdge->vexid_linked_to);
            curEdge=curEdge->next_edge;
        }
        while(!helperstack.empty()){
            helperval=helperstack.top();
            helperstack.pop();
            s.push(helperval);
        }
    }std::cout<<"\n";
}


/// calcu. without a stru.
int data_zr_hannoi_1[3] = { 2,0,0 };
void zr_hannoi(int num, int& start, int& tmp, int& end) {
	if (num == 1)
	{
		start--; end++;
		printf("%d  %d  %d\n", data_zr_hannoi_1[0], data_zr_hannoi_1[1], data_zr_hannoi_1[2]);
	}
	else {
		zr_hannoi(num - 1, start, end, tmp);
		zr_hannoi(1, start, tmp, end);
		zr_hannoi(num - 1, tmp, start, end);
	}
}


