#include "reflection.h"
#include "float.h"
#include "fstream"
#include "strstream"
#include <algorithm>
#define METHOD_3
map<int,int> face_patch;
temp *ref_table;
GBSolid gbsolid;
DiscretSolid discretsolid;
FergusonSurface *data_surface;
FergusonCurve *data_curve;
void datainitail(char *gm3file,HYBRID_MESH& mesh){
    for(int i=0;i<mesh.NumTris;++i){
        for(int j=0;j<3;++j)
        mesh.nodes[mesh.pTris[i].vertices[j]].flag=1;
    }
    read_gm3(gm3file,&gbsolid);
    vertex_G=gbsolid.vertex_G;
    GBCurve *curves=gbsolid.curve_G;
    GBFace *faces=gbsolid.face_G;
    int nf=gbsolid.nloops_G;
    int nc=gbsolid.ncurves_G;
    data_surface = new  FergusonSurface[nf];
    data_curve = new  FergusonCurve[nc];
    vertex *_vertex=new vertex[mesh.NumNodes];
    for(int i=0;i<nc;++i){
    createFergusonCurve(&gbsolid,&curves[i],&data_curve[i]);
    }
    for(int i=0;i<nf;++i){
    createFergusonFace(faces[i],&data_surface[i]);
    }
    double u,v;
    Vector coord;
    double distance,cls_u;
    double vol;
    int k;
    ref_table=new temp[mesh.NumNodes];
    double min=1e10;
    discretsolid.NumFacets=mesh.NumTris;
    discretsolid.NumPoints=mesh.NumNodes;
    discretsolid.discreteFacets=new DiscretFacet[discretsolid.NumFacets];
    discretsolid.discretPoints=new DiscretPoint[discretsolid.NumPoints];
    for(int i=0;i<discretsolid.NumFacets;++i){
        discretsolid.discretPoints[i].index=i;
               discretsolid.discretPoints[i].x=mesh.nodes[i].coord.x;
               discretsolid.discretPoints[i].y=mesh.nodes[i].coord.y;
               discretsolid.discretPoints[i].z=mesh.nodes[i].coord.z;
    }

    for(int i=0;i<discretsolid.NumFacets;++i){
        discretsolid.discreteFacets[i].index=i;
        discretsolid.discreteFacets[i].points[0]=mesh.pTris[i].vertices[0];
        discretsolid.discreteFacets[i].points[1]=mesh.pTris[i].vertices[1];
        discretsolid.discreteFacets[i].points[2]=mesh.pTris[i].vertices[2];
    }
    buildRelationshipByPoint(&discretsolid);
    buildFacetRelationshipByEdge(&discretsolid);



    for(int i=0;i<discretsolid.NumPoints;++i){
         if(discretsolid.discretPoints[i].linkedPoints.size()!=discretsolid.discretPoints[i].linkedFacets.size()){
             discretsolid.discretPoints[i].linkedPoints.erase(0);  //problem?
         }
    }

    Vector coord_1;
    Vector coord_2;
    for(int i=0;i<discretsolid.NumPoints;++i){
         coord_1=mesh.nodes[discretsolid.discretPoints[i].index].coord;
         for(auto j=discretsolid.discretPoints[i].linkedPoints.begin();j!=discretsolid.discretPoints[i].linkedPoints.end();++j)
         {coord_2=mesh.nodes[*j].coord;
         distance=coord_1.getDistance(coord_2);
         if(distance<min&&distance!=0)
         {
             min=distance;
         }
         }
    }
    vol=min/10;
#ifdef METHOD_1
              cout<<"method_1"<<endl;

           cout<<"start find curve!!"<<endl;
             for(int i=0;i<nc;++i){
                 for(int j=0;j<mesh.NumNodes;++j){
                     data_curve[i].project(mesh.nodes[j].coord,&coord,&cls_u);
                     distance=mesh.nodes[j].coord.getDistance(coord);
                     if(distance<vol){
                         if(_vertex[j].on_curve==false){
                             _vertex[j].curve=&data_curve[i];
                             _vertex[j].on_curve=true;
                             _vertex[j].flag=i;
                         }
                      }
                 }
             }
        cout<<"finished,start find face"<<endl;
             for(int i=0;i<nf;++i){
                 for(int j=0;j<mesh.NumNodes;++j){
                         {
                         data_surface[i].project(mesh.nodes[j].coord,&u,&v);
                         data_surface[i].param_to_coord(u,v,&coord);
                         if(coord.getDistance(mesh.nodes[j].coord)<vol&&_vertex[j].flag==-1)
                         {
                             _vertex[j].flag=i;
                             _vertex[j].surface=&data_surface[i];
                         }
                     }
                }
             }
             cout<<"finished"<<endl;
#endif

#ifdef METHOD_2
     cout<<"method_2"<<endl;
         double x[2],y[2],z[2];
         for(int i=0;i<nc;++i){
             if(data_curve[i].flag==0){
                 data_curve[i].get_box(x,y,z);
                 for(int j=0;j<mesh.NumNodes;++j){
                     if(mesh.nodes[j].flag==0)
                         continue;
                    if(mesh.nodes[j].coord.x>=x[0]-vol&&mesh.nodes[j].coord.y>=y[0]-vol&&mesh.nodes[j].coord.z>=z[0]-vol&&mesh.nodes[j].coord.x<=x[1]+vol&&mesh.nodes[j].coord.y<=y[1]+vol&&mesh.nodes[j].coord.z<=z[1]+vol){
                        data_curve[i].project(mesh.nodes[j].coord,&coord,&cls_u);
                        if(coord.getDistance(mesh.nodes[j].coord)<vol/10){
                            _vertex[j].on_curve=true;
                            _vertex[j].curve=&data_curve[i];
                            _vertex[j].flag=i;
                        }
                    }
                 }
             }
         }

         for(int i=0;i<nf;++i){
             if(data_surface[i].flag==0){
                 data_surface[i].get_box(x,y,z);
                 for(int j=0;j<mesh.NumNodes;++j){
                    if(_vertex[j].on_curve==false&&mesh.nodes[j].flag&&mesh.nodes[j].coord.x>=x[0]-vol&&mesh.nodes[j].coord.y>=y[0]-vol&&mesh.nodes[j].coord.z>=z[0]-vol&&mesh.nodes[j].coord.x<=x[1]+vol&&mesh.nodes[j].coord.y<=y[1]+vol&&mesh.nodes[j].coord.z<=z[1]+vol){
                        data_surface[i].project(mesh.nodes[j].coord,&u,&v);
                        data_surface[i].param_to_coord(u,v,&coord);
                        if(coord.getDistance(mesh.nodes[j].coord)<vol/100){
                            _vertex[j].surface = &data_surface[i];
                            _vertex[j].flag=i;
                        }
                    }
                 }
             }
         }
         cout<<"start violent exhaustion"<<endl;
        for(int i=0;i<nf;++i){
       //     cout<<i<<endl;
            for(int j=0;j<mesh.NumNodes;++j){
                    if(_vertex[j].flag==-1){
                    data_surface[i].project(mesh.nodes[j].coord,&u,&v);
                    data_surface[i].param_to_coord(u,v,&coord);
                    if(coord.getDistance(mesh.nodes[j].coord)<vol)
                    {
                        _vertex[j].flag=i;
                        _vertex[j].surface=&data_surface[i];
                    }
                }
              }
        }
#endif



#ifdef METHOD_3
         cout<<"method_3"<<endl;
         cout<<"start Start building relationships of discrete and continuous"<<endl;;
        vector<vector<int>> curve_head;
        vector<int> temp;
        bool headfind=0,endfind=0;
        for(int i=0;i<nc;++i){
            for(int j=0;j<mesh.NumNodes;++j){
                if(mesh.nodes[j].flag==0)
                    continue;
                if(data_curve[i].start_coord().getDistance(mesh.nodes[j].coord)<vol)
                {
                    headfind=1;
                    temp.push_back(j);
                }
                if(data_curve[i].end_coord().getDistance(mesh.nodes[j].coord)<vol){
                    temp.push_back(j);
                    endfind=1;
                }
            }
            if(headfind&&endfind){
                if(_vertex[temp[0]].on_curve==false){
                _vertex[temp[0]].on_curve=true;
                _vertex[temp[0]].flag=i;
                _vertex[temp[0]].curve=&data_curve[i];
                }
                if(_vertex[temp[1]].on_curve==false){
                _vertex[temp[1]].on_curve=true;
                _vertex[temp[1]].flag=i;
                _vertex[temp[1]].curve=&data_curve[i];
                }
                data_curve[i].flag=1;
            }
            else
  //              cout<<"curve "<<i<<" find error"<<endl;
                headfind=0;
                endfind=0;
            curve_head.push_back(temp);
            temp.clear();
        }

//        for(int i=0;i<curve_head.size();++i){
//            cout<<curve_head[i].size()<<endl;
//        }

        vector<vector<int>> curve_node;
        int start_node;
        int end_node;
        int pre_node=-1;
          k=0;
        DiscretPoint current_node;
        for(int i=0;i<nc;++i){
        temp.clear();
        start_node=curve_head[i][0];
        end_node=curve_head[i][1];
        current_node=discretsolid.discretPoints[start_node];
        temp.push_back(start_node);
        while(current_node.index!=end_node){
           for(auto j=current_node.linkedPoints.begin();j!=current_node.linkedPoints.end();++j){
               data_curve[i].project(mesh.nodes[*j].coord,&coord,&cls_u);
               distance=mesh.nodes[*j].coord.getDistance(coord);
               if((*j)!=pre_node&&distance<vol/*prolbem*/){
               temp.push_back(*j);
               _vertex[*j].on_curve=true;
               _vertex[*j].curve=&data_curve[i];
               _vertex[*j].flag=i;
               pre_node=current_node.index;
               current_node=discretsolid.discretPoints[*j];
               break;
           }
        }
        }
        pre_node=-1;
        curve_node.push_back(temp);
        }
//       cout<<"curve_node"<<curve_node.size()<<endl;
        vector<set<int>> iner_loop;
        vector<set<int>> loop;
        set<int> _temp;
        int flag=1;
        for(int i=0;i<gbsolid.nloops_G;++i){
            _temp.clear();
            for(int j=0;j<gbsolid.loop_G[i].nlc;++j){
                if(data_curve[gbsolid.loop_G[i].cp_t[j]->index].flag==0)
               {
                    flag=0;
     //           cout<<i<<"  is not loop!"<<endl;
                }
            }
            if(flag){
             for(int j=0;j<gbsolid.loop_G[i].nlc;++j){
                for(int k=0;k<curve_node[gbsolid.loop_G[i].cp_t[j]->index].size();++k){
                    int n=curve_node[gbsolid.loop_G[i].cp_t[j]->index][k];
                    _temp.insert(n);
                }
            }
             data_surface[i].flag=1;
         }
            flag=1;
            loop.push_back(_temp);
      }
    //    cout<<"loop :"<<loop.size()<<endl;
        for(int i=0;i<gbsolid.nloops_G;++i){
            _temp.clear();
            for(auto j=loop[i].begin();j!=loop[i].end();++j) {
                current_node=discretsolid.discretPoints[*j];
               discretsolid.discretPoints[*j].curvePoint=true;
                 for(auto k=current_node.linkedPoints.begin();k!=current_node.linkedPoints.end();++k)
                  {
                        if(loop[i].find(*k)==loop[i].end()){
                     data_surface[i].project(mesh.nodes[*k].coord,&u,&v);
                    data_surface[i].param_to_coord(u,v,&coord);
                    distance=mesh.nodes[*k].coord.getDistance(coord);
                    if(distance<vol)
                       {
                        discretsolid.discretPoints[*k].curvePoint=true;
                        _vertex[*k].surface=&data_surface[i];
                        _vertex[*k].flag=i;
                        _temp.insert(*k);
                       }
                  }
                 }
            }
            iner_loop.push_back(_temp);
        }
 //       cout<<"iner_loop "<<iner_loop.size()<<endl;


        vector<vector<int>> face_node;
        vector<int> face_node_temp;
        for(int i=0;i<gbsolid.nloops_G;++i){
            for(auto j=iner_loop[i].begin();j!=iner_loop[i].end();++j){
                current_node=discretsolid.discretPoints[*j];
                for(auto k=current_node.linkedPoints.begin();k!=current_node.linkedPoints.end();++k){
                    if(!discretsolid.discretPoints[*k].curvePoint){
                        current_node=discretsolid.discretPoints[*k];
                    //    face_node_temp.push_back(*k);
                        break;
                    }
                }
           while(!discretsolid.discretPoints[current_node.index].curvePoint){
               face_node_temp.push_back(current_node.index);
               _vertex[current_node.index].surface=&data_surface[i];
               _vertex[current_node.index].flag=i;
               discretsolid.discretPoints[current_node.index].curvePoint=true;
               for(auto k=current_node.linkedPoints.begin();k!=current_node.linkedPoints.end();++k){
                   if(!discretsolid.discretPoints[*k].curvePoint){
                       current_node=discretsolid.discretPoints[*k];
                       break;
                   }
               }
           }
            }
           face_node.push_back(face_node_temp);
           face_node_temp.clear();
        }
         cout<<"start violent exhaustion"<<endl;
        for(int i=0;i<nf;++i){
       //     cout<<i<<endl;
            for(int j=0;j<mesh.NumNodes;++j){
                    if(_vertex[j].flag==-1){
                    data_surface[i].project(mesh.nodes[j].coord,&u,&v);
                    data_surface[i].param_to_coord(u,v,&coord);
                    if(coord.getDistance(mesh.nodes[j].coord)<vol)
                    {
                        _vertex[j].flag=i;
                        _vertex[j].surface=&data_surface[i];
                    }
                }
              }
        }
#endif

 #ifdef METHOD_4
     cout<<"method_4"<<endl;
     cout<<"start find End point of curve"<<endl;
    vector<vector<int>> curve_head;
    vector<int> temp;
    bool headfind=0,endfind=0;
    for(int i=0;i<nc;++i){
        for(int j=0;j<mesh.NumNodes;++j){
            if(mesh.nodes[j].flag==0)
                continue;
            if(data_curve[i].start_coord().getDistance(mesh.nodes[j].coord)<vol){
                headfind=1;
                temp.push_back(j);
            }
            if(data_curve[i].end_coord().getDistance(mesh.nodes[j].coord)<vol){
                temp.push_back(j);
                endfind=1;
            }
        }
        if(headfind&&endfind){
            _vertex[temp[0]].on_curve=true;
            _vertex[temp[0]].flag=i;
            _vertex[temp[0]].curve=&data_curve[i];
            _vertex[temp[1]].on_curve=true;
            _vertex[temp[1]].flag=i;
            _vertex[temp[1]].curve=&data_curve[i];
            data_curve[i].flag=1;
        }
            headfind=0;
            endfind=0;
        curve_head.push_back(temp);
        temp.clear();
    }

    for(int i=0;i<curve_head.size();++i){
    }
    cout<<endl;
     cout<<"use box find remain curve"<<endl;
    vector<vector<int>> box_nodes;
    double x[2],y[2],z[2];
    for(int i=0;i<nc;++i){
        temp.clear();
        if(data_curve[i].flag==0){
            data_curve[i].get_box(x,y,z);
            for(int j=0;j<mesh.NumNodes;++j){
                if(mesh.nodes[j].flag==0)
                    continue;
               if(mesh.nodes[j].coord.x>=x[0]-vol&&mesh.nodes[j].coord.y>=y[0]-vol&&mesh.nodes[j].coord.z>=z[0]-vol&&mesh.nodes[j].coord.x<=x[1]+vol&&mesh.nodes[j].coord.y<=y[1]+vol&&mesh.nodes[j].coord.z<=z[1]+vol){
                   data_curve[i].project(mesh.nodes[j].coord,&coord,&cls_u);
                   if(coord.getDistance(mesh.nodes[j].coord)<vol){
                       temp.push_back(j);
                       _vertex[j].on_curve=true;
                       _vertex[j].curve=&data_curve[i];
                       _vertex[j].flag=i;
                   }
               }
            }
            box_nodes.push_back(temp);
        }
    }
    cout<<endl;
    cout<<"use End point to find curve"<<endl;
    vector<vector<int>> curve_node;
    int start_node;
    int end_node;
    int pre_node=-1;
      k=0;
    DiscretPoint current_node;
    for(int i=0;i<nc;++i){

    if(!data_curve[i].flag)
      {
        if(!box_nodes[k].empty())
          //  data_curve[i].flag=1;     //problem
        curve_node.push_back(box_nodes[k++]);
        continue;
    }
    temp.clear();
    start_node=curve_head[i][0];
    end_node=curve_head[i][1];
    current_node=discretsolid.discretPoints[start_node];
    temp.push_back(start_node);
   // cout<<"start_node"<<start_node<<endl;
    //cout<<mesh.nodes[start_node]
    while(current_node.index!=end_node){
       for(auto j=current_node.linkedPoints.begin();j!=current_node.linkedPoints.end();++j){
           data_curve[i].project(mesh.nodes[*j].coord,&coord,&cls_u);
           distance=mesh.nodes[*j].coord.getDistance(coord);
           if((*j)!=pre_node&&distance<vol/*prolbem*/){
           temp.push_back(*j);
           _vertex[*j].on_curve=true;
           _vertex[*j].curve=&data_curve[i];
           _vertex[*j].flag=i;
           pre_node=current_node.index;
           current_node=discretsolid.discretPoints[*j];
           break;
       }
    }
    }

    pre_node=-1;
    curve_node.push_back(temp);
    }
    cout<<"find loop "<<endl;
    vector<set<int>> iner_loop;
    vector<set<int>> loop;
    set<int> _temp;
    int flag=1;
    for(int i=0;i<gbsolid.nloops_G;++i){
        _temp.clear();
        for(int j=0;j<gbsolid.loop_G[i].nlc;++j){
            if(data_curve[gbsolid.loop_G[i].cp_t[j]->index].flag==0)
           {
                flag=0;
            cout<<i<<"  is not loop!"<<endl;
            }
        }
        if(flag){
         for(int j=0;j<gbsolid.loop_G[i].nlc;++j){
            for(int k=0;k<curve_node[gbsolid.loop_G[i].cp_t[j]->index].size();++k){
                int n=curve_node[gbsolid.loop_G[i].cp_t[j]->index][k];
                _temp.insert(n);
            }
        }
         data_surface[i].flag=1;
     }
        flag=1;
        loop.push_back(_temp);
  }
    cout<<"loop :"<<loop.size()<<endl;
    for(int i=0;i<gbsolid.nloops_G;++i){
        cout<<i<<endl;
        _temp.clear();
        for(auto j=loop[i].begin();j!=loop[i].end();++j) {
            current_node=discretsolid.discretPoints[*j];
           discretsolid.discretPoints[*j].curvePoint=true;
             for(auto k=current_node.linkedPoints.begin();k!=current_node.linkedPoints.end();++k)
              {
                    if(loop[i].find(*k)==loop[i].end()){
                 data_surface[i].project(mesh.nodes[*k].coord,&u,&v);
                data_surface[i].param_to_coord(u,v,&coord);
                distance=mesh.nodes[*k].coord.getDistance(coord);
                if(distance<vol)
                   {
                    discretsolid.discretPoints[*k].curvePoint=true;
                    _vertex[*k].surface=&data_surface[i];
                    _vertex[*k].flag=i;
                    _temp.insert(*k);
                   }
              }
             }
        }
        iner_loop.push_back(_temp);
    }
    cout<<"iner_loop "<<iner_loop.size()<<endl;
    box_nodes.clear();
    cout<<"use box find remain face"<<endl;
    for(int i=0;i<nf;++i){
        temp.clear();
        if(data_surface[i].flag==0){
            data_surface[i].get_box(x,y,z);
            for(int j=0;j<mesh.NumNodes;++j){
               if(mesh.nodes[j].flag&&mesh.nodes[j].coord.x>=x[0]-vol&&mesh.nodes[j].coord.y>=y[0]-vol&&mesh.nodes[j].coord.z>=z[0]-vol&&mesh.nodes[j].coord.x<=x[1]+vol&&mesh.nodes[j].coord.y<=y[1]+vol&&mesh.nodes[j].coord.z<=z[1]+vol){
                   data_curve[i].project(mesh.nodes[j].coord,&coord,&cls_u);
                   if(_vertex[j].on_curve==false&&coord.getDistance(mesh.nodes[j].coord)<vol){
                       temp.push_back(j);
                       _vertex[j].surface = &data_surface[i];
                       _vertex[j].flag=i;
                   }
               }
            }
            box_nodes.push_back(temp);
        }
    }

    vector<vector<int>> face_node;
    vector<int> face_node_temp;
    for(int i=0;i<gbsolid.nloops_G;++i){
        cout<<i<<endl;
        if(!data_curve[i].flag)
          {
            if(!box_nodes[k].empty())
                data_curve[i].flag=1;
            curve_node.push_back(box_nodes[k++]);
            continue;
        }
        for(auto j=iner_loop[i].begin();j!=iner_loop[i].end();++j){
            current_node=discretsolid.discretPoints[*j];
            for(auto k=current_node.linkedPoints.begin();k!=current_node.linkedPoints.end();++k){
                if(!discretsolid.discretPoints[*k].curvePoint){
                    current_node=discretsolid.discretPoints[*k];
                    goto flag;
                }
            }

     }
        flag:
       while(!discretsolid.discretPoints[current_node.index].curvePoint){
           face_node_temp.push_back(current_node.index);
           _vertex[current_node.index].surface=&data_surface[i];
           _vertex[current_node.index].flag=i;
           discretsolid.discretPoints[current_node.index].curvePoint=true;
           for(auto k=current_node.linkedPoints.begin();k!=current_node.linkedPoints.end();++k){
               if(!discretsolid.discretPoints[*k].curvePoint){
                   current_node=discretsolid.discretPoints[*k];
                //   face_node_temp.push_back(*k);
                   break;
               }
           }
       }
       face_node.push_back(face_node_temp);
       face_node_temp.clear();
    }
#endif

    for(int i=0;i<mesh.NumNodes;++i){
        if(_vertex[i].on_curve==true){
             ref_table[i].curve_id=_vertex[i].flag;
             ref_table[i].patch_id=-1;
             ref_table[i].node_id=i;
        }
        else {
            ref_table[i].patch_id=_vertex[i].flag;
            ref_table[i].curve_id=-1;
            ref_table[i].node_id=i;
       }
    }


    cout<<"reflection success!"<<endl;
    ofstream file;
    file.open("reflection_table.txt");
    if(file.is_open()){
        file<<"node_id"<<"  "<<"curve_id"<<"  "<<"patch_id"<<endl;
        for(int i=0;i<mesh.NumNodes;++i){
            file<<ref_table[i].node_id<<"  \t\t"<<ref_table[i].curve_id<<"  \t\t"<<ref_table[i].patch_id<<endl;
        }
        file.close();
    }
}


Refletion::Refletion()
{

}
void Refletion::initial(char *gm3file,HYBRID_MESH& mesh)
{
    datainitail(gm3file,mesh);
    for(int i=0;i<mesh.NumTris;++i){
        int patch_id1=ref_table[mesh.pTris[i].vertices[0]].patch_id;
        int patch_id2=ref_table[mesh.pTris[i].vertices[1]].patch_id;
        int patch_id3=ref_table[mesh.pTris[i].vertices[2]].patch_id;
        if(patch_id1!=-1&&patch_id2!=-1&&patch_id3!=-1)
        {
            face_patch[i]=patch_id1;
        }
        else if(patch_id1==-1&&patch_id2==-1&&patch_id3==-1)
     {
            cout<<"tris "<<i<<endl;
            cout<<"  "<<mesh.pTris[i].vertices[0]<<"   "<<mesh.pTris[i].vertices[1]<<"   "<<mesh.pTris[i].vertices[2]<<endl;
          Vector vec=(mesh.nodes[mesh.pTris[i].vertices[0]].coord+mesh.nodes[mesh.pTris[i].vertices[1]].coord+mesh.nodes[mesh.pTris[i].vertices[2]].coord)/3;
          face_patch[i]=findsurface(ref_table[mesh.pTris[i].vertices[0]].curve_id,ref_table[mesh.pTris[i].vertices[1]].curve_id,ref_table[mesh.pTris[i].vertices[2]].curve_id,vec);
     }
        else if(patch_id1==-1){
            if(patch_id2==-1)
                face_patch[i]=patch_id3;
            else if(patch_id3==-1)
                face_patch[i]=patch_id2;
            else
                face_patch[i]=patch_id2;
        }
        else if(patch_id2==-1){
            if(patch_id1==-1)
                face_patch[i]=patch_id3;
            else if(patch_id3==-1)
                face_patch[i]=patch_id1;
            else
            face_patch[i]=patch_id1;
        }
        else if(patch_id3==-1){
            if(patch_id1==-1)
                face_patch[i]=patch_id2;
            else if(patch_id2==-1)
                face_patch[i]=patch_id1;
            else
            face_patch[i]=patch_id1;
        }

    }

    for(int i=0;i<mesh.NumTris;++i){
        if(face_patch[i]!=-1)
      {  gbsolid.face_G[face_patch[i]].facets.insert(i);
        mesh.pTris[i].iSurface=face_patch[i];
        }
    }
    this->curves=data_curve;
    this->sufaces=data_surface;
    set<int>::iterator it;
    for(int i=0;i<gbsolid.nloops_G;i++){
        it=gbsolid.face_G[i].facets.begin();
        for(int j=0;j<gbsolid.face_G[i].facets.size();++j){

            if(it!=gbsolid.face_G[i].facets.end())
           {
            if(subject_table[*(it)]==0)
            subject_table[*(it)]=gbsolid.face_G[i].index+1;//problem
            it++;
           }
        }
        }

}

Vector Refletion::subject_test(int patch_ID_1, int patch_ID_2, Vector coord)
{
    double u,v;
    Vector n=coord;
    this->sufaces[subject_table[patch_ID_1]].project(coord,&u,&v);
    return n;
}

void Refletion::attach_face(int face_ID, int patch_ID)
{
    subject_table[face_ID]=patch_ID;
}

int Refletion::detach_face(int face_ID)
{
     int result =(*subject_table.find(face_ID)).second;
     subject_table.erase(face_ID);
     return result;
}
 int count=0;
 Vector Refletion::subject_patch_id(int patch_ID_1,int patch_ID_2,Vector coord){
     Vector cls_pnt=coord;
     double u,v,cls_u;
     if(patch_ID_1!=patch_ID_2)
     {
//         cout<<"patcd_id_1 "<<patch_ID_1<<" patch_ID_2 "<<patch_ID_2<<endl;
          FergusonCurve *temp;
        temp=findcurve(patch_ID_2,patch_ID_1,coord);
               if(temp==nullptr)
               {
  //                 cout<<patch_ID_1<<"   "<<patch_ID_2<<endl;
                   cout<<"can't find curve!!!!!!!"<<endl;
               }
           else
               temp->project(coord,&cls_pnt,&cls_u);//problem
     }
     else
     {
         this->sufaces[patch_ID_1].project(coord,&u,&v);
         this->sufaces[patch_ID_1].param_to_coord(u,v,&cls_pnt);
     }
     return cls_pnt;
 }

Vector Refletion::subject_face_id(int face_ID_1,int face_ID_2,Vector coord)
{
   Vector cls_pnt=coord;
   double u,v,cls_u;
   if(subject_table[face_ID_1]!=subject_table[face_ID_2])
   {
    // cout<<"face_ID_1: "<<face_ID_1<<"face_ID_1: "<<face_ID_2<<endl;
        FergusonCurve *temp;
      temp=findcurve(subject_table[face_ID_1]-1,subject_table[face_ID_2]-1,coord);
             if(temp==nullptr)
             {
                 cout<<"can't find curve!!!!!!!"<<endl;
                 cout<<"face_ID_1: "<<face_ID_1<<"face_ID_1: "<<face_ID_2<<endl;
                 cout<<"patch_id_1="<<subject_table[face_ID_1]<<"patch_id_2="<<subject_table[face_ID_2]<<endl;
             }
         else
              temp->project(coord,&cls_pnt,&cls_u);//problem
   }
   else
   {
       this->sufaces[subject_table[face_ID_1]-1].project(coord,&u,&v);
       this->sufaces[subject_table[face_ID_1]-1].param_to_coord(u,v,&cls_pnt);
   }
   return cls_pnt;
}


FergusonCurve* Refletion::findcurve(int patch_id_1, int patch_id_2,Vector coord)
{
    double min=1e10;
    double distance=0;
    FergusonCurve *savecurve=new FergusonCurve;
    savecurve=nullptr;
    FergusonCurve *temp=new FergusonCurve;
    for(int i=0;i<gbsolid.loop_G[patch_id_1].nlc;++i){
       for(int j=0;j<gbsolid.loop_G[patch_id_2].nlc;++j){
    if(gbsolid.loop_G[patch_id_1].cp_t[i]->index==gbsolid.loop_G[patch_id_2].cp_t[j]->index)
    {
        temp=&data_curve[gbsolid.loop_G[patch_id_1].cp_t[i]->index];
        Vector cls_pnt;
        double cls_u;
        temp->project(coord,&cls_pnt,&cls_u);
        distance=cls_pnt.getDistance(coord);
        if( min>distance ){
            min=distance;
            savecurve=temp;
           }

         }
       }
    }
    return savecurve;
}

int Refletion::findsurface(int curve_id_1, int curve_id_2, int curve_id_3,Vector vector)
{
    int a=0,b=0,c=0;
    double u,v;
    Vector coord;
    for(int i=0;i<gbsolid.nloops_G;++i){
        for(int j=0;j<gbsolid.loop_G[i].nlc;++j)
        {
            if(gbsolid.loop_G[i].cp_t[j]->index==curve_id_1)
                             a=1;
            if(gbsolid.loop_G[i].cp_t[j]->index==curve_id_2)
                             b=1;
            if(gbsolid.loop_G[i].cp_t[j]->index==curve_id_3)
                             c=1;
        }
        if(a&&b&&c)
        {
            data_surface[i].project(vector,&u,&v);
            data_surface[i].param_to_coord(u,v,&coord);
            if(coord.getDistance(vector)<0.01)
               {
                cout<<"surface : "<<i<<endl;
                return i;
            }
        }
       else  if(i==gbsolid.nloops_G-1)
        {
            cout<<"1   can't find surface!!!!!!!"<<endl;
            return -1;

        }
        a=0;b=0;c=0;
    }
    cout<<"2   can't find surface!!!!!!!"<<endl;
    return -1;
}

void Refletion::edge_to_face(HYBRID_MESH &mesh,int &face_id_1, int &face_id_2,string edge_name)
{
        vector<string>::iterator strIter;
        vector<int> face_id;
        strIter=lines.begin();
           while(strIter!=lines.end()){
               strIter=find(strIter,lines.end(),edge_name);
               if(strIter!=lines.end()){
                   face_id.push_back((strIter-lines.begin())/3);//problem?strIter
                   strIter++;
               }


           }//#add

    if(face_id.size()>1){
    for(int k=0;k<face_id.size()-1;++k){
        if(mesh.pTris[face_id[k]].iSurface!=mesh.pTris[face_id[k+1]].iSurface){
            face_id_1=face_id[k];
            face_id_2=face_id[k+1];
            break;
        }
        face_id_1=face_id[k];
        face_id_2=face_id[k];
    }
    }    //simple judge
 //#add
}

void Refletion::initial_edge_table(HYBRID_MESH &mesh)
{
    int v[3];
    lines.clear();
    for(int i=0;i<mesh.NumTris;i++)
    {
        v[0]=mesh.pTris[i].vertices[0];
        v[1]=mesh.pTris[i].vertices[1];
        v[2]=mesh.pTris[i].vertices[2];
        sort(v,v+3);
        lines.push_back(IntToString(v[0])+"_"+IntToString(v[1]));
        lines.push_back(IntToString(v[1])+"_"+IntToString(v[2]));
        lines.push_back(IntToString(v[0])+"_"+IntToString(v[2]));
    }
}

