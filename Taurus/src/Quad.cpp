#include "Quad.h"
void Quad::add(float x, float y, int idx) {
    size++;
    float tempx = mass[0], tempy = mass[1];
    bool sameflag = false;
    float mult1 = (float)(size - 1.0f) / size;
    float mult2 = (float)1.0f / size;
    mass[0] = mass[0] * mult1 + mult2 * x;
    mass[1] = mass[1] * mult1 + mult2 * y;
    if (is_leaf & !have_node) {
        have_node = true;
        nodeindex = idx;
        return;
    }

    if (is_leaf) {
        if ((x - mass[0]) * (x - mass[0]) + (y - mass[1]) * (y - mass[1]) > 1e-12) {
            children[0] = new Quad(center[0] - 0.25f * w, center[1] - 0.25f * h, 0.5f * w, 0.5f * h);
            children[1] = new Quad(center[0] + 0.25f * w, center[1] - 0.25f * h, 0.5f * w, 0.5f * h);
            children[2] = new Quad(center[0] - 0.25f * w, center[1] + 0.25f * h, 0.5f * w, 0.5f * h);
            children[3] = new Quad(center[0] + 0.25f * w, center[1] + 0.25f * h, 0.5f * w, 0.5f * h);
            is_leaf = false;
            if (mult1 <= 0.0f) {
                printf("In here must more than one node and mult1>0.0f");
                exit(1);
            }
            int index = 0;
            if (tempx > center[0])index += 1;
            if (tempy > center[1])index += 2;
            children[index]->add(tempx, tempy, nodeindex);

        }
        else {
            children[0] = new Quad(center[0] - 0.25f * w, center[1] - 0.25f * h, 0.5f * w, 0.5f * h);
            children[1] = new Quad(center[0] + 0.25f * w, center[1] - 0.25f * h, 0.5f * w, 0.5f * h);
            children[2] = new Quad(center[0] - 0.25f * w, center[1] + 0.25f * h, 0.5f * w, 0.5f * h);
            children[3] = new Quad(center[0] + 0.25f * w, center[1] + 0.25f * h, 0.5f * w, 0.5f * h);
            is_leaf = false;
            if (mult1 <= 0.0f) {
                printf("In here must more than one node and mult1>0.0f");
                exit(1);
            }
            int index = 0;
            if (tempx > center[0])index += 1;
            if (tempy > center[1])index += 2;
            children[index]->add(tempx, tempy, nodeindex);
            sameflag = true;
        }
    }

    if (!is_leaf) {
        int index = 0;
        if (x > center[0])index += 1;
        if (y > center[1])index += 2;
        index = (index + sameflag) % 4;
        children[index]->add(x, y, idx);
    }
}
Quad::Quad(std::vector<Node>& nodes, int N) {
    float maxx = -1e9, minx = 1e9, maxy = -1e9, miny = 1e9;
    for (int i = 0; i < N; i++) {
        maxx = nodes[i].x > maxx ? nodes[i].x : maxx;
        minx = nodes[i].x < minx ? nodes[i].x : minx;
        maxy = nodes[i].y > maxy ? nodes[i].y : maxy;
        miny = nodes[i].y < miny ? nodes[i].y : miny;
    }
    is_leaf = true;
    have_node = false;
    nodeindex = -1;
    size = 0;
    mass[0] = mass[1] = .0;
    center[0] = (maxx + minx) / 2.0f;
    center[1] = (maxy + miny) / 2.0f;

    w = (maxx - minx) > (maxy - miny) ? (maxx - minx) : (maxy - miny);
    w = 1.01f * w;
    h = w;
    children = (Quad**)malloc(4 * sizeof(Quad*));
    children[0] = children[1] = children[2] = children[3] = NULL;
    for (int i = 0; i < N; i++) {
        add(nodes[i].x, nodes[i].y, i);
    }

}

Quad::Quad(float cx, float cy, float width, float height) {
    is_leaf = true;
    have_node = false;
    size = 0;
    mass[0] = mass[1] = .0;
    center[0] = cx, center[1] = cy;
    w = width;
    h = height;
    nodeindex = -1;
    children = (Quad**)malloc(4 * sizeof(Quad*));
    children[0] = children[1] = children[2] = children[3] = NULL;
}

Quad::~Quad() {
    for (int i = 0; i < 4; i++) {
        if (children[i] != NULL)
            delete children[i];
    }
    free(children);
}

void QuadForce(Quad* qt, float& fx, float& fy, float posx, float posy, float theta, float alpha, int nodeindex,float dp,int layer) {
    //fprintf(stderr,"Hello QuadForce\n");
    //bool debugFlag = false;
    float eps=0;
    if(dp>1.0){eps=0.1;}else{eps=0.01;}
    if (qt->size <= 0.0)return;
    float mvx = (posx - qt->mass[0]);
    float mvy = (posy - qt->mass[1]);
    float dist = mvx * mvx + mvy * mvy;
    // if(debugFlag)fprintf(stderr,"%.16f %f %f\n",dist,mvx,mvy);
    if (qt->is_leaf && dist <= 1e-4) {
        if (qt->size > 1) {
            printf("No qtnode have more than one node");
            exit(1);
        }
        //if (debugFlag)fprintf(stderr, "%d %d\n", qt->nodeindex, nodeindex);
        if (qt->nodeindex != nodeindex) {
            if (mvx == 0.0) { mvx = 0.001f * (rand() / RAND_MAX - 0.5f), dist += mvx * mvx; }
            if (mvy == 0.0) { mvy = 0.001f * (rand() / RAND_MAX - 0.5f), dist += mvy * mvy; }
            if (dist < 1e-12) {
                dist = sqrtf(0.01f * dist);
            }
            // if(debugFlag)fprintf(stderr,"dist,%f %f %f %f %f %f %f %f %f\n",qt->mass[0],qt->mass[1],qt->w,qt->size,dist,mvx,mvy,fx,fy);
            //********* force define  *********//
            dist=sqrtf(dist) + eps;
            fx -= alpha * mvx * (qt->size - 1) / pow(dist,dp+1);
            fy -= alpha * mvy * (qt->size - 1) / pow(dist,dp+1);
        }
    }
    else if (qt->is_leaf || dist / (qt->w * qt->w) > theta * theta) {

        if (mvx == 0.0) { mvx = 0.001f * (rand() / RAND_MAX - 0.5f), dist += mvx * mvx; }
        if (mvy == 0.0) { mvy = 0.001f * (rand() / RAND_MAX - 0.5f), dist += mvy * mvy; }
        if (dist <= 1e-5) {
            dist = sqrtf(0.01f * dist) + eps;
        }
        // if(debugFlag)fprintf(stderr,"theta,%f %f %f %f %f %f %f %f %f\n",qt->mass[0],qt->mass[1],qt->w,qt->size,dist,mvx,mvy,fx,fy);
        //********* force define  *********//
        dist=sqrtf(dist) + eps;
        fx -= alpha * mvx * qt->size /pow(dist,dp+1);
        fy -= alpha * mvy * qt->size /pow(dist,dp+1);
    }
    else {
        layer++;
       if(layer>12){return;}
        QuadForce(qt->children[0], fx, fy, posx, posy, theta, alpha, nodeindex,dp,layer);
        QuadForce(qt->children[1], fx, fy, posx, posy, theta, alpha, nodeindex,dp,layer);
        QuadForce(qt->children[2], fx, fy, posx, posy, theta, alpha, nodeindex,dp,layer);
        QuadForce(qt->children[3], fx, fy, posx, posy, theta, alpha, nodeindex,dp,layer);
    }
}
void Quad_t_Force(Quad* qt, float& fx, float& fy, float posx, float posy, float theta, float alpha, int nodeindex,float dp,float tf,int layer) {
    //fprintf(stderr,"Hello QuadForce\n");
    //bool debugFlag = false;
    float eps=0;
    //if(dp>1.0){eps=0.1;}else{eps=0.01;}
    if (qt->size <= 0.0)return;
    float mvx = (posx - qt->mass[0]);
    float mvy = (posy - qt->mass[1]);
    float dist = mvx * mvx + mvy * mvy;
    //cout<<"dp="<<dp<<" tf="<<tf<<endl;
    // if(debugFlag)fprintf(stderr,"%.16f %f %f\n",dist,mvx,mvy);
    if (qt->is_leaf && dist <= 1e-4) {
        if (qt->size > 1) {
            printf("No qtnode have more than one node");
            exit(1);
        }
        //if (debugFlag)fprintf(stderr, "%d %d\n", qt->nodeindex, nodeindex);
        if (qt->nodeindex != nodeindex) {
            if (mvx == 0.0) { mvx = 0.001f * (rand() / RAND_MAX - 0.5f), dist += mvx * mvx; }
            if (mvy == 0.0) { mvy = 0.001f * (rand() / RAND_MAX - 0.5f), dist += mvy * mvy; }
            if (dist < 1e-12) {
                dist = sqrtf(0.01f * dist);
            }
            // if(debugFlag)fprintf(stderr,"dist,%f %f %f %f %f %f %f %f %f\n",qt->mass[0],qt->mass[1],qt->w,qt->size,dist,mvx,mvy,fx,fy);
            //********* force define  *********//
            dist=sqrtf(dist) + eps;
            float t_dist=1.0/(1+dist*dist);
            fx -= alpha * mvx * (qt->size - 1)*pow(t_dist,tf) *pow(dist,dp)/dist;
            fy -= alpha * mvy * (qt->size - 1)*pow(t_dist,tf) *pow(dist,dp)/dist;
        }
    }
    else if (qt->is_leaf || dist / (qt->w * qt->w) > theta * theta) {

        if (mvx == 0.0) { mvx = 0.001f * (rand() / RAND_MAX - 0.5f), dist += mvx * mvx; }
        if (mvy == 0.0) { mvy = 0.001f * (rand() / RAND_MAX - 0.5f), dist += mvy * mvy; }
        if (dist <= 1e-5) {
            dist = sqrtf(0.01f * dist) + eps;
        }
        // if(debugFlag)fprintf(stderr,"theta,%f %f %f %f %f %f %f %f %f\n",qt->mass[0],qt->mass[1],qt->w,qt->size,dist,mvx,mvy,fx,fy);
        //********* force define  *********//
        dist=sqrtf(dist) + eps;
        float t_dist=1.0/(1+dist*dist);
        fx -= alpha * mvx * qt->size*pow(t_dist,tf) *pow(dist,dp)/dist;
        fy -= alpha * mvy * qt->size*pow(t_dist,tf) *pow(dist,dp)/dist;
    }
    else {
        layer++;
       if(layer>8){return;}
        Quad_t_Force(qt->children[0], fx, fy, posx, posy, theta, alpha, nodeindex,dp,tf,layer);
        Quad_t_Force(qt->children[1], fx, fy, posx, posy, theta, alpha, nodeindex,dp,tf,layer);
        Quad_t_Force(qt->children[2], fx, fy, posx, posy, theta, alpha, nodeindex,dp,tf,layer);
        Quad_t_Force(qt->children[3], fx, fy, posx, posy, theta, alpha, nodeindex,dp,tf,layer);
    }
}