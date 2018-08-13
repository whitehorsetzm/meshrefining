#ifndef TABLE_H
#define TABLE_H
#include<map>
#include<vector>
#include<set>
#include<queue>
#include"dataclass.h"
#include"dataio.h"
#include"ferguson_curve.h"
#include"ferguson_surface.h"
#include"adaptive_io.h"
using namespace std;
using namespace EEMAS;
extern FergusonSurface *data_surface;
extern FergusonCurve *data_curve;
extern GBSolid gbsolid;
extern DiscretSolid dssolid;
extern HYBRID_MESH mesh;
void datainitail(char *gm3file,HYBRID_MESH& mesh);
class vertex{
public:
    vertex(){
        neighbor.clear();
        flag=-1;
        surface=nullptr;
        curve=nullptr;
        index=-1;
        face.clear();
        on_curve=false;
    }
  int index;
  set<int> neighbor;
  set<int> face;  //face contain this pointer
  FergusonSurface *surface;
  FergusonCurve   *curve;
  int flag;
  bool on_curve;
};
class table
{
public:
    table();
    table(table *a){
        this->subject_table=a->subject_table;
        this->curves=a->curves;
        this->sufaces=a->sufaces;
    }
    void initial(char *gm3file,HYBRID_MESH& mesh);//��ĳ�ʼ������ͨ��subject()������ȡһ������ÿ��������ľ��루���ߵ���һ��������Ļ��ľ��룩��ѡȡ����ġ��������Ҫͨ��������Ĺ�ϵ�õ�
    //������Ĺ�ϵ������������Ҫһ�����ݽṹ�������������ȴ�����
    void attach_face(int face_ID,int patch_ID);
    int detach_face(int face_ID);
    Vector subject(int face_ID_1,int face_ID_2,Vector coord);
    FergusonCurve* findcurve(int patch_id_1,int patch_id_2,Vector coord);
     map<int,int> subject_table;
    FergusonSurface *sufaces;
    FergusonCurve   *curves;
};
extern table _table;
class temp{
public:
    int patch_id;
    int curve_id;
    int node_id;
    temp():patch_id(-1),curve_id(-1),node_id(-1){/*cout<<"sucess"<<endl;*/}
};
//int vn=mesh.NumNodes;
//vertex* vertices=new vertex[vn];//�±�Ϊindex
//void init(){
//    for(int i=0;i<elemsize;i++){
//        for(int j=0;j<3;++j){
//            vertices[mesh.pTris[i].vertices[j]].index=mesh.pTris[i].vertices[j];
//            vertices[mesh.pTris[i].vertices[j]].neighbor.insert(mesh.pTris[i].vertices[(j+1)%3]);
//            vertices[mesh.pTris[i].vertices[j]].neighbor.insert(mesh.pTris[i].vertices[(j+2)%3]);
//            vertices[mesh.pTris[i].vertices[j]].face.insert(mesh.pTris[i].index);
//        }
//    }
//}

//class loop_{
//  set<int>            node;//index
//  FergusonCurve     *curve;//loop�ɼ���curve���
//};
//class face_{
//  set<int>            node;
//  FergusonSurface    *suface;
//};
//loop_ edges[nloop];
//for(int i=0;i<n_node;++i){
//    node[i]
//}
//ÿ��loop_�ﴢ���ż����������߹��ɵĻ�����Щ���������е������.
//������뼸�����ߵĹ�ϵ֮�󣬿��Եõ������һȦ���㣬������Ȧ��������ɫ�㷨���Եõ�һ���������Ӧ�����������档
//��ɫ�㷨��Ҫ���ݽṹ��֧��
//�����ݽṹ��Ҫ�����ڽӵ����Ϣ������ڽ���Ϣ����ͨ����������������õ���

//vertex vetices[vn] �±�Ϊ�����ȫ��INDEX������������Ƭ����3��������в�����������neighbor�м������������㣬��neighbor��û���������㣬��ֱ�Ӽ��룬����������
//����ͼ�㷨��������ʱ������ǰ����Ϊ�߽�㣬����Ҫ�ж���һ�������Ƿ�����Ҫ�ļ������ϣ�����ǣ����ǲ����ʣ�������ǣ�����������
//����߽�㣻
//��μ���ĳ�����Ƿ���ĳ�������棿 ���õ��򼸺���ͶӰ��Ȼ��ͶӰ�����붥�����������롣�����ԣ�

//face_ face*=new face_[nf];
//for(int i=0;i<nloops;++i){
//    queue<int> q;
//    q.push(edges[i].node[0]);
//    while(!q.empty()){
//        int n=q.front();
//        face[i].node.insert(n);
//        q.pop();
//        vertices[n].flag=true;
//    for(int i=0;i<vertices[n].neighbor.size();++i){
//        if(!edges[i].node.find(t))
//        {
//            //����ҵ���˵���ڱ߽���
//            //��Ҫ�ж��¸����������ڻ�������
//            if(//����&&vertices[n].flag==false)
//                    )
//               q.push(vertices[n].neighbor[i]);
//        else{
//            if(vertices[n].flag==false)
//            q.push(vertices[n].neighbor[i]);
//        }
//    }
//    }

//}
// }
//   face�����б�����ÿ������������Ӧ�������
//   edges�����б�����ÿ�������߶�Ӧ�������
//��ȫ��������֮��Ҫת��Ϊ������Ĺ�ϵ�����ͨ�������棿
//���´�mesh�������濪ʼ������ʱ�����е������ĸ�������ߵ���Ϣ���������������ͺܿ��ˡ�
//�������˳���ȶ�mesh��������б������ж�ÿ������3�����Ƿ��е�����ĳ��


#endif // TABLE_H
