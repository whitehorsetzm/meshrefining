#include "reflection.h"
#include "float.h"
#include "fstream"
#include "strstream"
//table _table;
map<int,int> face_patch;
temp *ref_table;
GBSolid gbsolid;
DiscretSolid discretsolid;
FergusonSurface *data_surface;
FergusonCurve *data_curve;
//void getmaxcoord_curve(FergusonCurve curve,double &x,double &y,double &z){
//    curve.
//}
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
    \
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
    int n=-1,m=-1;
    for(int i=0;i<discretsolid.NumPoints;++i){
         coord_1=mesh.nodes[discretsolid.discretPoints[i].index].coord;
         for(auto j=discretsolid.discretPoints[i].linkedPoints.begin();j!=discretsolid.discretPoints[i].linkedPoints.end();++j)
         {coord_2=mesh.nodes[*j].coord;
         distance=coord_1.getDistance(coord_2);
         if(distance<min)
         {
             min=distance;
             n=i;
             m=*j;
         }
         }
    }
    cout<<n<<"  "<<m<<endl;
    vol=min/10;
     cout<<"vol ="<<vol<<endl;
    //suan fa 3
//    vector<vector<int>> curve_head;
//    vector<int> temp;
//    bool headfind=0,endfind=0;
//    for(int i=0;i<nc;++i){
//        for(int j=0;j<mesh.NumNodes;++j){
//            if(mesh.nodes[j].flag==0)
//                continue;
//            if(data_curve[i].start_coord().getDistance(mesh.nodes[j].coord)<vol){
//                headfind=1;
//                temp.push_back(j);
//            }
//            if(data_curve[i].end_coord().getDistance(mesh.nodes[j].coord)<vol){
//                temp.push_back(j);
//                endfind=1;
//            }
//        }
////        if(!headfind)
////        cout<<data_curve[i].start_coord().x<<" "<<data_curve[i].start_coord().y<<" "<<data_curve[i].start_coord().z<<endl;
////        if(!endfind)
////        cout<<data_curve[i].end_coord().x<<" "<<data_curve[i].start_coord().y<<" "<<data_curve[i].start_coord().z<<endl;
//        if(headfind&&endfind){
//            _vertex[temp[0]].on_curve=true;
//            _vertex[temp[0]].flag=i;
//            _vertex[temp[0]].curve=&data_curve[i];
//            _vertex[temp[1]].on_curve=true;
//            _vertex[temp[1]].flag=i;
//            _vertex[temp[1]].curve=&data_curve[i];
//            data_curve[i].flag=1;
//       //     cout<<"find curve!!!!!!!"<<endl;
//        }
//            headfind=0;
//            endfind=0;
//        curve_head.push_back(temp);
//        temp.clear();
//    }

//    for(int i=0;i<curve_head.size();++i){
//        cout<<curve_head[i].size()<<endl;
//    }
//    vector<vector<int>> box_nodes;
//    double x[2],y[2],z[2];
//    for(int i=0;i<nc;++i){
//        temp.clear();
//        if(data_curve[i].flag==0){
//            data_curve[i].get_box(x,y,z);
////            cout<<"curve :"<<i<<endl;
////            cout<<"min x"<<x[0]<<"  "<<x[1]<<endl;
////            cout<<"min y"<<y[0]<<"  "<<y[1]<<endl;
////            cout<<"min z"<<z[0]<<"  "<<z[1]<<endl;
//            for(int j=0;j<mesh.NumNodes;++j){
//                if(mesh.nodes[j].flag==0)
//                    continue;
//               if(mesh.nodes[j].coord.x>=x[0]-vol&&mesh.nodes[j].coord.y>=y[0]-vol&&mesh.nodes[j].coord.z>=z[0]-vol&&mesh.nodes[j].coord.x<=x[1]+vol&&mesh.nodes[j].coord.y<=y[1]+vol&&mesh.nodes[j].coord.z<=z[1]+vol){
//                   data_curve[i].project(mesh.nodes[j].coord,&coord,&cls_u);
//                 //  cout<<j<<endl;
//                   if(coord.getDistance(mesh.nodes[j].coord)<vol){
//                       cout<<"add a pointer"<<endl;
//                       temp.push_back(j);
//                       _vertex[j].on_curve=true;
//                       _vertex[j].curve=&data_curve[i];
//                       _vertex[j].flag=i;
//                   }
//               }
//            }
//            box_nodes.push_back(temp);
//        }
//    }
//    for(int i=0;i<nc;++i){
//        if(!data_curve[i].flag)
//        cout<<i<<endl;
//    }
//    cout<<"!!!!!!!!!!!!!!!!!!!!!!!"<<box_nodes.size()<<endl;
//    vector<vector<int>> curve_node;
//    int start_node;
//    int end_node;
//    int pre_node=-1;
//      k=0;
//    DiscretPoint current_node;
//    for(int i=0;i<nc;++i){

//    if(!data_curve[i].flag)
//      {
//        if(!box_nodes[k].empty())
//          //  data_curve[i].flag=1;     //problem
//        curve_node.push_back(box_nodes[k++]);
//        continue;
//    }
//         cout<<"curve "<<i<<endl;
//    temp.clear();
//    start_node=curve_head[i][0];
//    end_node=curve_head[i][1];
//    current_node=discretsolid.discretPoints[start_node];
//    temp.push_back(start_node);
//   // cout<<"start_node"<<start_node<<endl;
//    //cout<<mesh.nodes[start_node]
//    while(current_node.index!=end_node){
//       for(auto j=current_node.linkedPoints.begin();j!=current_node.linkedPoints.end();++j){
//   //        cout<<"*"<<endl;
//           data_curve[i].project(mesh.nodes[*j].coord,&coord,&cls_u);
//           distance=mesh.nodes[*j].coord.getDistance(coord);
//       //    cout<<pre_node<<"   "<<current_node.index<<endl;
//           if((*j)!=pre_node&&distance<vol/*prolbem*/){
//           temp.push_back(*j);
//           _vertex[*j].on_curve=true;
//           _vertex[*j].curve=&data_curve[i];
//           _vertex[*j].flag=i;
//           pre_node=current_node.index;
//           current_node=discretsolid.discretPoints[*j];
//           break;
//       }
//    }
//    }

//    pre_node=-1;
//    curve_node.push_back(temp);
//  //  cout<<i<<endl;
//    }
//   cout<<"curve_node"<<curve_node.size()<<endl;

//    vector<set<int>> iner_loop;
//    vector<set<int>> loop;
//    set<int> _temp;
//    int flag=1;
//    for(int i=0;i<gbsolid.nloops_G;++i){
//        _temp.clear();
//        for(int j=0;j<gbsolid.loop_G[i].nlc;++j){
//            if(data_curve[gbsolid.loop_G[i].cp_t[j]->index].flag==0)
//           {
//                flag=0;
//            cout<<i<<"  is not loop!"<<endl;
//            }
//        }
//        if(flag){
//         for(int j=0;j<gbsolid.loop_G[i].nlc;++j){
//            for(int k=0;k<curve_node[gbsolid.loop_G[i].cp_t[j]->index].size();++k){
//                int n=curve_node[gbsolid.loop_G[i].cp_t[j]->index][k];
//                _temp.insert(n);
//            }
//        }
//         data_surface[i].flag=1;
//     }
//        flag=1;
//        loop.push_back(_temp);
//  }
//    cout<<"loop :"<<loop.size()<<endl;
//    for(int i=0;i<gbsolid.nloops_G;++i){
//        cout<<i<<endl;
//        _temp.clear();
//        for(auto j=loop[i].begin();j!=loop[i].end();++j) {
//            current_node=discretsolid.discretPoints[*j];
//           discretsolid.discretPoints[*j].curvePoint=true;
//             for(auto k=current_node.linkedPoints.begin();k!=current_node.linkedPoints.end();++k)
//              {
//                    if(loop[i].find(*k)==loop[i].end()){
//                 data_surface[i].project(mesh.nodes[*k].coord,&u,&v);
//                data_surface[i].param_to_coord(u,v,&coord);
//                distance=mesh.nodes[*k].coord.getDistance(coord);
//                if(distance<vol)
//                   {
//                    discretsolid.discretPoints[*k].curvePoint=true;
//                    _vertex[*k].surface=&data_surface[i];
//                    _vertex[*k].flag=i;
//                    _temp.insert(*k);
//                   }
//              }
//             }
//        }
//        iner_loop.push_back(_temp);
//    }
//    cout<<"iner_loop "<<iner_loop.size()<<endl;
//    box_nodes.clear();
//    for(int i=0;i<nf;++i){
//        temp.clear();
//        if(data_surface[i].flag==0){
//            data_surface[i].get_box(x,y,z);
//            cout<<"surface :"<<i<<endl;
////            cout<<"min x"<<x[0]<<"  "<<x[1]<<endl;
////            cout<<"min y"<<y[0]<<"  "<<y[1]<<endl;
////            cout<<"min z"<<z[0]<<"  "<<z[1]<<endl;
//            for(int j=0;j<mesh.NumNodes;++j){
//               if(mesh.nodes[j].flag&&mesh.nodes[j].coord.x>=x[0]-vol&&mesh.nodes[j].coord.y>=y[0]-vol&&mesh.nodes[j].coord.z>=z[0]-vol&&mesh.nodes[j].coord.x<=x[1]+vol&&mesh.nodes[j].coord.y<=y[1]+vol&&mesh.nodes[j].coord.z<=z[1]+vol){
//                   data_curve[i].project(mesh.nodes[j].coord,&coord,&cls_u);
//                 //  cout<<j<<endl;
//                   if(_vertex[j].on_curve==false&&coord.getDistance(mesh.nodes[j].coord)<vol){
//           //            cout<<"add a pointer"<<endl;
//                       temp.push_back(j);
//                       _vertex[j].surface = &data_surface[i];
//                       _vertex[j].flag=i;
//                   }
//               }
//            }
//            box_nodes.push_back(temp);
//        }
//    }


//    vector<vector<int>> face_node;
//    vector<int> face_node_temp;
//    for(int i=0;i<gbsolid.nloops_G;++i){
//        if(!data_curve[i].flag)
//          {
//            if(!box_nodes[k].empty())
//                data_curve[i].flag=1;
//            curve_node.push_back(box_nodes[k++]);
//            continue;
//        }
//        for(auto j=iner_loop[i].begin();j!=iner_loop[i].end();++j){
//            current_node=discretsolid.discretPoints[*j];
//            for(auto k=current_node.linkedPoints.begin();k!=current_node.linkedPoints.end();++k){
//                if(!discretsolid.discretPoints[*k].curvePoint){
//                    current_node=discretsolid.discretPoints[*k];
//                //    face_node_temp.push_back(*k);
//               //     cout<<"first"<<*k<<endl;
//               //     cout<<current_node.index<<endl;
//                    goto flag;
//                }
//            }

//     }
//        flag:
//       while(!discretsolid.discretPoints[current_node.index].curvePoint){
//    //       cout<<current_node.index<<endl;
//           face_node_temp.push_back(current_node.index);
//           _vertex[current_node.index].surface=&data_surface[i];
//           _vertex[current_node.index].flag=i;
//           discretsolid.discretPoints[current_node.index].curvePoint=true;
//           for(auto k=current_node.linkedPoints.begin();k!=current_node.linkedPoints.end();++k){
//               if(!discretsolid.discretPoints[*k].curvePoint){
//                   current_node=discretsolid.discretPoints[*k];
//                //   face_node_temp.push_back(*k);
//                   break;
//               }
//           }
//       }
//       face_node.push_back(face_node_temp);
//       face_node_temp.clear();
//    }
////    for(int i=0;i<gbsolid.nloops_G;++i){
////        for(int c=0;c<gbsolid.loop_G[i].nlc;++c){
////            for(int j=0;j<curve_node[c].size();++j){
////                ref_table[curve_node[c][j]].curve_id=c;
////                 ref_table[curve_node[c][j]].patch_id=-1;
////                  ref_table[curve_node[c][j]].node_id=curve_node[c][j];
////            }
////        }
////        for(auto j=iner_loop[i].begin();j!=iner_loop[i].end();++j){
////            ref_table[*j].curve_id=-1;
////             ref_table[*j].patch_id=i;
////              ref_table[*j].node_id=*j;
////        }
////        for(int j=0;j<face_node[i].size();++j){
////            ref_table[face_node[i][j]].curve_id=-1;
////             ref_table[face_node[i][j]].patch_id=i;
////              ref_table[face_node[i][j]].node_id=face_node[i][j];
////        }
////    }
////    for(int i=0;i<mesh.NumNodes;++i){
////    //    cout<<_vertex[i].flag<<endl;
////    }

         double x[2],y[2],z[2];
         for(int i=0;i<nc;++i){
             if(data_curve[i].flag==0){
                 data_curve[i].get_box(x,y,z);
     //            cout<<"curve :"<<i<<endl;
     //            cout<<"min x"<<x[0]<<"  "<<x[1]<<endl;
     //            cout<<"min y"<<y[0]<<"  "<<y[1]<<endl;
     //            cout<<"min z"<<z[0]<<"  "<<z[1]<<endl;
                 for(int j=0;j<mesh.NumNodes;++j){
                     if(mesh.nodes[j].flag==0)
                         continue;
                    if(mesh.nodes[j].coord.x>=x[0]-vol&&mesh.nodes[j].coord.y>=y[0]-vol&&mesh.nodes[j].coord.z>=z[0]-vol&&mesh.nodes[j].coord.x<=x[1]+vol&&mesh.nodes[j].coord.y<=y[1]+vol&&mesh.nodes[j].coord.z<=z[1]+vol){
                        data_curve[i].project(mesh.nodes[j].coord,&coord,&cls_u);
                      //  cout<<j<<endl;
                        if(coord.getDistance(mesh.nodes[j].coord)<vol/10){
                           // cout<<"add a pointer"<<endl;
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
                 cout<<"surface :"<<i<<endl;
                 cout<<"min x"<<x[0]<<"  "<<x[1]<<endl;
                 cout<<"min y"<<y[0]<<"  "<<y[1]<<endl;
                 cout<<"min z"<<z[0]<<"  "<<z[1]<<endl;
                 for(int j=0;j<mesh.NumNodes;++j){
                    if(_vertex[j].on_curve==false&&mesh.nodes[j].flag&&mesh.nodes[j].coord.x>=x[0]-vol&&mesh.nodes[j].coord.y>=y[0]-vol&&mesh.nodes[j].coord.z>=z[0]-vol&&mesh.nodes[j].coord.x<=x[1]+vol&&mesh.nodes[j].coord.y<=y[1]+vol&&mesh.nodes[j].coord.z<=z[1]+vol){
                        data_surface[i].project(mesh.nodes[j].coord,&u,&v);
                        data_surface[i].param_to_coord(u,v,&coord);
                      //    cout<<j<<endl;
                        if(coord.getDistance(mesh.nodes[j].coord)<vol/100){
                       //     cout<<"add a pointer"<<endl;
                            _vertex[j].surface = &data_surface[i];
                            _vertex[j].flag=i;
                        }
                    }
                 }
             }
         }


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
//    // suanfa 1
//    vertex *_vertex=new vertex[mesh.NumNodes];
//    for(int i=0;i<mesh.NumTris;++i){
//        for(int j=0;j<3;++j)
//        {
//            //cout<<mesh.pTris[i].vertices[j]<<endl;
//            _vertex[mesh.pTris[i].vertices[j]].neighbor.insert(mesh.pTris[i].vertices[(j+1)%3]);
//            _vertex[mesh.pTris[i].vertices[j]].neighbor.insert(mesh.pTris[i].vertices[(j+2)%3]);

//        }
//    }
//    for(int i=0;i<mesh.NumNodes;++i){
//        for(auto a=_vertex[i].neighbor.begin();a!=_vertex[i].neighbor.end();++a){
//      //      cout<<"!!!!"<<*a<<endl;
//        }
//    }

//    vector<vector<int>> curve_node;
//    for(int i=0;i<nc;++i){
//        vector<int> temp;
//        for(int j=0;j<mesh.NumNodes;++j){
//            data_curve[i].project(mesh.nodes[j].coord,&coord,&cls_u);
//            distance=mesh.nodes[j].coord.getDistance(coord);
//            if(distance<vol){
//                temp.push_back(j);
//                if(_vertex[j].on_curve==false){
//                    _vertex[j].curve=&data_curve[i];
//                    _vertex[j].on_curve=true;
//                    _vertex[j].flag=i;
//                }
//             }
//        }
//        curve_node.push_back(temp);
//  //      cout<<"i="<<i<<"   "<<temp.size()<<endl;
//    }

//    for(int i=0;i<curve_node.size();++i){
//        for(int j=0;j<curve_node[i].size();++j){
//      //           cout<<"curve_ID= "<<i<<"node_ID= "<<curve_node[i][j]<<endl;
//        }
//    }
//    for(int i=0;i<curve_node.size();++i){
//    for(int j=0;j<curve_node[i].size();++j){
//    //      cout<<curve_node[i][j]<<"  ";
//    }
//   // cout<<endl;
//    }
//     vector<set<int>> loop_node;
//    for(int i=0;i<nf;++i){
//        set<int> temp;
//        for(int j=0;j<gbsolid.loop_G[i].nlc;++j){
//      //      cout<<"the "<<gbsolid.loop_G[i].cp_t[j]->index<<" curve "<<"has "<<curve_node[gbsolid.loop_G[i].cp_t[j]->index].size()<<endl;
//            for(int k=0;k<curve_node[gbsolid.loop_G[i].cp_t[j]->index].size();++k)
//            {
//                temp.insert(curve_node[gbsolid.loop_G[i].cp_t[j]->index][k]);
//              //    cout<<curve_node[gbsolid.loop_G[i].cp_t[j]->index][k]<<endl;
//            }

//        }
//        loop_node.push_back(temp);
//    }
//        for(int i=0;i<loop_node.size();++i){
//            for(auto j=loop_node[i].begin();j!=loop_node[i].end();++j){
//       //              cout<<"face_ID= "<<i<<"node_ID= "<<*j<<endl;
//            }
//        }


//      //  for(int i=0;i<mesh.NumNodes;++i)
//      //           cout<<"_vertex[i]"<<_vertex[i].neighbor.size()<<endl;
////for(int i=0;i<loop_node.size();++i){
////    for(auto j=loop_node[i].begin();j!=loop_node[i].end();++j){
////        for(auto k=_vertex[*j].neighbor.begin();k!=_vertex[*j].neighbor.end();++k)
////            {
//// //           cout<<*j<<endl;
//////           cout<<"k= "<<*k<<endl;
////            data_surface[i].project(mesh.nodes[*k].coord,&u,&v);
////            data_surface[i].param_to_coord(u,v,&coord);
//////                cout<<"i"<<i<<endl;
//////            if(i=-52)
//////                cout<<coord.getDistance(mesh.nodes[*k].coord)<<"?????????????????"<<endl;    //problem;
////            if(coord.getDistance(mesh.nodes[*k].coord)<vol&&_vertex[*k].flag==-1)
////            {
////                _vertex[*k].flag=i;
////                _vertex[*k].surface=&data_surface[i];
////            }
//////            cout<<"test";
////        }
////   }
////}
////for(int i=0;i<loop_node.size();++i){
////    for(auto j=loop_node[i].begin();j!=loop_node[i].end();++j){
////        for(auto k=_vertex[*j].neighbor.begin();k!=_vertex[*j].neighbor.end();++k)
////            {
//// //           cout<<*j<<endl;
//////           cout<<"k= "<<*k<<endl;
////            data_surface[i].project(mesh.nodes[*k].coord,&u,&v);
////            data_surface[i].param_to_coord(u,v,&coord);
//////                cout<<"i"<<i<<endl;
//////            if(i=-52)
//////                cout<<coord.getDistance(mesh.nodes[*k].coord)<<"?????????????????"<<endl;    //problem;
////            if(coord.getDistance(mesh.nodes[*k].coord)<vol&&_vertex[*k].flag==-1)
////            {
////                _vertex[*k].flag=i;
////                _vertex[*k].surface=&data_surface[i];
////            }
//////            cout<<"test";
////        }
////   }
////}
////vertex *current;

////for(int i=0;i<mesh.NumNodes;++i){
////    int t=-1;
////    if(_vertex[i].flag==-1)
////     current=&_vertex[i];
////    while(current!=nullptr&&current->flag==-1){
////        if(!current->neighbor.empty()){
////        for(auto k=current->neighbor.begin();k!=current->neighbor.end();++k){
////            if(t==-1&&_vertex[*k].flag!=-1)
////                t=_vertex[*k].flag;
////        }

////        }
////        if(t==-1)
////            break;
////        current->flag=t;
////        for(auto k=current->neighbor.begin();k!=current->neighbor.end();++k){
////            if(_vertex[*k].flag==-1)
////              {
////                current=&_vertex[*k];
////                break;
////            }
////        }
////    }
////    t=-1;
////    cout<<i<<endl;
////}

//for(int i=0;i<nf;++i){
//    for(int j=0;j<mesh.NumNodes;++j){
//            {
// //           cout<<*j<<endl;
////           cout<<"k= "<<*k<<endl;
//            data_surface[i].project(mesh.nodes[j].coord,&u,&v);
//            data_surface[i].param_to_coord(u,v,&coord);
////                cout<<"i"<<i<<endl;
//            if(coord.getDistance(mesh.nodes[j].coord)<vol&&_vertex[j].flag==-1)
//            {
//                _vertex[j].flag=i;
//                _vertex[j].surface=&data_surface[i];
//            }
////            cout<<"test";
//        }
//   }
//}

////cout<<"test"<<endl;
//for(int i=0;i<mesh.NumNodes;++i){
////    cout<<_vertex[i].flag<<endl;
//}

//for(int i=0;i<mesh.NumNodes;++i){
//    if(_vertex[i].on_curve==true){
//         ref_table[i].curve_id=_vertex[i].flag;
//         ref_table[i].patch_id=-1;
//         ref_table[i].node_id=i;
//    }
//    else {
//        ref_table[i].patch_id=_vertex[i].flag;
//        ref_table[i].curve_id=-1;
//        ref_table[i].node_id=i;
//   }

//}





//      //suanfa 2
//    for(int i=0;i<mesh.NumNodes;++i){
//        min==DBL_MAX;
//        k=0;
//        for(int j=0;j<nf;++j){
//            data_surface[j].project(mesh.nodes[i].coord,&u,&v);
//            data_surface[j].param_to_coord(u,v,&coord);
//            distance=mesh.nodes[i].coord.getDistance(coord);
//            if(distance<=vol){
//                f[m]=j;
//                m++;
//                k=j;
//            }
//           if(distance<min){
//                min=distance;
//                k=j;
//            }
//        }
//      //  cout<<"distance:"<<min<<endl;
//      //   cout<<m<<endl; //distance has problem
//        if(m>=2){
//            min=DBL_MAX;
//            for(int j=0;j<nc;++j){
//                data_curve[j].project(mesh.nodes[i].coord,&coord,&cls_u);
//                distance=mesh.nodes[i].coord.getDistance(coord);
//                if(distance<min){
//                     min=distance;
//                     k=j;
//                 }
//            }
//            data_curve[k].project(mesh.nodes[i].coord,&coord,&cls_u);
//            distance=mesh.nodes[i].coord.getDistance(coord);
////             cout<<"k="<<k<<"curve:"<<distance<<endl;

//            ref_table[i].curve_id=k;
//            ref_table[i].node_id=i;
//        }
//        else{

//            data_surface[k].project(mesh.nodes[i].coord,&u,&v);
//            data_surface[k].param_to_coord(u,v,&coord);
//            distance=mesh.nodes[i].coord.getDistance(coord);
//          //  cout<<"k="<<k<<"face:"<<distance<<endl;

//             ref_table[i].patch_id=k;
//             ref_table[i].node_id=i;
//        }
//         m=0;
//         min=DBL_MAX;
//    }
//for(int i=0;i<nc;++i){
//    data_curve[i]
//}






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
}  //#ok


table::table()
{

}
void table::initial(char *gm3file,HYBRID_MESH& mesh)
{
    datainitail(gm3file,mesh);
 //   ref_table=new temp[mesh.NumNodes];
//    while(file>>s)){
//        for(int i=0;i<3;++i)
//        {
//        }
//    }
    for(int i=0;i<mesh.NumTris;++i){
        if(mesh.nodes[i].flag==0)
            continue;
        int patch_id1=ref_table[mesh.pTris[i].vertices[0]].patch_id;
        int patch_id2=ref_table[mesh.pTris[i].vertices[1]].patch_id;
        int patch_id3=ref_table[mesh.pTris[i].vertices[2]].patch_id;
        //cout<<patch_id1<<" "<<patch_id2<<" "<<patch_id3<<endl;
        if(patch_id1!=-1&&patch_id2!=-1&&patch_id3!=-1)
        {
            face_patch[i]=patch_id1;
        }
        else if(patch_id1==-1&&patch_id2==-1&&patch_id3==-1)
     {
          Vector vec=(mesh.nodes[mesh.pTris[i].vertices[0]].coord+mesh.nodes[mesh.pTris[i].vertices[1]].coord+mesh.nodes[mesh.pTris[i].vertices[2]].coord)/3;
          face_patch[i]=findsurface(ref_table[mesh.pTris[i].vertices[0]].curve_id,ref_table[mesh.pTris[i].vertices[1]].curve_id,ref_table[mesh.pTris[i].vertices[2]].curve_id,vec);

         cout<<ref_table[mesh.pTris[i].vertices[0]].curve_id<<"  "<<ref_table[mesh.pTris[i].vertices[1]].curve_id<<" "<<ref_table[mesh.pTris[i].vertices[2]].curve_id<<endl;
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

           // cout<<"patch_id "<<gbsolid.face_G[i].index<<"face id "<<*(it)<<endl;
            if(subject_table[*(it)]==0)
            subject_table[*(it)]=gbsolid.face_G[i].index+1;//problem
            it++;
           }
        }
        }

}

Vector table::subject_test(int patch_ID_1, int patch_ID_2, Vector coord)
{
    double u,v;
    Vector n=coord;
    this->sufaces[subject_table[patch_ID_1]].project(coord,&u,&v);
    return n;
}

void table::attach_face(int face_ID, int patch_ID)
{
    subject_table[face_ID]=patch_ID;
}

int table::detach_face(int face_ID)
{
     int result =(*subject_table.find(face_ID)).second;
     subject_table.erase(face_ID);
     return result;
}
 int count=0;
 Vector table::subject_patch_id(int patch_ID_1,int patch_ID_2,Vector coord){
     Vector cls_pnt=coord;
     double u,v,cls_u;
     if(patch_ID_1!=patch_ID_2)
     {
          FergusonCurve *temp;
        temp=findcurve(patch_ID_2,patch_ID_1,coord);
               if(temp==nullptr)
               {
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

Vector table::subject_face_id(int face_ID_1,int face_ID_2,Vector coord)
{
   Vector cls_pnt=coord;
   double u,v,cls_u;
   if(subject_table[face_ID_1]!=subject_table[face_ID_2])
   {
     //cout<<"face_ID_1: "<<face_ID_1<<"face_ID_1: "<<face_ID_2<<endl;
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

     //  cout<<"test";

   }
   else
   {
     //  cout<<"face_ID_1: "<<face_ID_1<<"face_ID_1: "<<face_ID_2<<endl;
     //  cout<<"patch_id_1="<<subject_table[face_ID_1]<<"patch_id_2="<<subject_table[face_ID_2]<<endl;
       this->sufaces[subject_table[face_ID_1]-1].project(coord,&u,&v);
       this->sufaces[subject_table[face_ID_1]-1].param_to_coord(u,v,&cls_pnt);
   }
   return cls_pnt;
}


FergusonCurve* table::findcurve(int patch_id_1, int patch_id_2,Vector coord)
{
    FergusonCurve *temp=new FergusonCurve;
  //  cout<<"patch_id_1="<<patch_id_1<<"patch_id_2="<<patch_id_2<<endl;
  //  cout<<gbsolid.loop_G[patch_id_1].nlc<<"  "<<gbsolid.loop_G[patch_id_2].nlc<<endl;
    for(int i=0;i<gbsolid.loop_G[patch_id_1].nlc;++i){
       for(int j=0;j<gbsolid.loop_G[patch_id_2].nlc;++j){
  //  cout<<"gbsolid.loop_G[patch_id_1].cp_t[i]="<<gbsolid.loop_G[patch_id_1].cp_t[i]->index<<"  gbsolid.loop_G[patch_id_2].cp_t[j]"<<gbsolid.loop_G[patch_id_2].cp_t[j]->index<<endl;
    if(gbsolid.loop_G[patch_id_1].cp_t[i]->index==gbsolid.loop_G[patch_id_2].cp_t[j]->index)
    {
        createFergusonCurve(&gbsolid,gbsolid.loop_G[patch_id_1].cp_t[i],temp);

        Vector cls_pnt;
        double cls_u;
        temp->project(coord,&cls_pnt,&cls_u);
     //   cout<<"distance = "<<cls_pnt.getDistance(coord)<<endl;
        if(cls_pnt.getDistance(coord)<0.05)         //problem;
               return temp;
    }
       }
    }
    return nullptr;
}

int table::findsurface(int curve_id_1, int curve_id_2, int curve_id_3,Vector vector)
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

