#include "stdDefinitions.h"
#include <math.h>

#include "agpWiegelmannPar.h"
#include "agmScalarField.h"
#include "agmVectorField.h"
#include "mfoGlobals.h"

#include "debug_data_trace_win.h"

#define DEBUG_OUT_PATH "g:\\temp\\"
class out_field
{
protected:
    int nProc;
    int iterN;
    int depth;
    char fn[1024], buff[32];

public:

    out_field(int _nProc, int _iterN, int _depth)
        : nProc(_nProc)
        , iterN(_iterN)
        , depth(_depth)
    {
    }

    void out(char *filename, CagmVectorField *v)
    {
#ifdef DEBUG_DATA_TRACE
        if (depth == 1)
        {
            strcpy(fn, DEBUG_OUT_PATH);
            strcat(fn, filename);
            strcat(fn, "_");
            itoa(iterN, buff, 32);
            strcat(fn, buff);
            strcat(fn, ".bin");
            debug_data_trace_win *d = new debug_data_trace_win(fn);
            d->write_vector(v);
            delete d;
        }
#endif
    }

    void out(char *filename, CagmScalarField *v)
    {
#ifdef DEBUG_DATA_TRACE
        if (depth == 1)
        {
            strcpy(fn, DEBUG_OUT_PATH);
            strcat(fn, filename);
            strcat(fn, "_");
            itoa(iterN, buff, 32);
            strcat(fn, buff);
            strcat(fn, ".bin");
            debug_data_trace_win *d = new debug_data_trace_win(fn);
            d->write_scalar(v);
            delete d;
        }
#endif
    }
};

//-----------------------------------------------------------------------
CagpWiegelmann::CagpWiegelmann(int *_N, int _n_threads, int _stencil
    , CagmVectorField *_sourceB, CagmScalarField *_sourceW
    , double *_vcos
    , CagmVectorFieldOps *_outF
    , CagmVectorField *_baseField, CagmVectorField *_baseWeight
    , CagmScalarField *_absField, CagmScalarField *_absWeight
    , CagmScalarField *_losField, CagmScalarField *_losWeight
    , CagmScalarField *_bottomWeight
    , int _priority
    )
    : CagmMetrics(1, _N[2])
    , vB(_sourceB)
    , sW(_sourceW)
    , vdircos(_vcos)
    , vF(_outF)
    , baseField(_baseField)
    , baseWeight(_baseWeight)
    , absField(_absField)
    , absWeight(_absWeight)
    , losField(_losField)
    , losWeight(_losWeight)
    , bottomWeight(_bottomWeight)
    , Bmod(nullptr)
    , Jmod(nullptr)
    , Binv(nullptr)
    , Jinv(nullptr)
    , JxB(nullptr)
    , st(nullptr)
    , priority(_priority)
{
    memcpy(N, vB->N, 3*sizeof(int));
    double *steps = vB->GetSteps();

    vgradW = new CagmVectorField(N);
    vgradW->SetSteps(steps);

    vgradW->grad(sW);

    int Nb[3];
    Nb[0] = N[0];
    Nb[1] = N[1];
    Nb[2] = 1;
    vBottom = new CagmVectorField(Nb);

    B2 = new CagmScalarField(N);
    B2->SetSteps(steps);
    rotB = new CagmVectorField(N);
    rotB->SetSteps(steps);
    divB = new CagmScalarField(N);
    divB->SetSteps(steps);
    Wa = new CagmVectorField(N);
    Wa->SetSteps(steps);
    Wb = new CagmVectorField(N);
    Wb->SetSteps(steps);
    Wa2Wb2 = new CagmScalarField(N);
    Wa2Wb2->SetSteps(steps);
    WaxB = new CagmVectorField(N);
    WaxB->SetSteps(steps);
    WbxB = new CagmScalarField(N);
    WbxB->SetSteps(steps);
    v = new CagmVectorField(N);
    v->SetSteps(steps);
    s = new CagmScalarField(N);
    s->SetSteps(steps);

    if (WiegelmannGetMetricsTheta)
    {
        st  =  new CagmScalarField(N);
        Binv = new CagmScalarField(N);
        Jinv = new CagmScalarField(N);
        Bmod = new CagmScalarField(N);
        Jmod = new CagmScalarField(N);
        JxB  = new CagmScalarField(N);
    }

    main_proc = new TaskQueueProcessor(_n_threads);
    nProc = main_proc->get_num_proc();

    factory = new PCOTaskFactory();
    supervisor = new PCOSupervisor(N, factory, main_proc->get_sync());

    for (int i = 0; i < nProc; i++)
        processors.push_back(new PCOProcessor(i, main_proc->get_sync(), this, _stencil));
}

//-----------------------------------------------------------------------
CagpWiegelmann::~CagpWiegelmann()
{
    for (int i = 0; i < nProc; i++)
        delete processors[i];

    delete supervisor;
    delete factory;
    delete main_proc;

    delete vgradW;

    delete vBottom;

    delete B2;
    delete rotB;
    delete divB;
    delete Wa;
    delete Wb;
    delete WaxB;
    delete WbxB;
    delete v;
    delete s;

    delete Binv;
    delete Jinv;
    delete Bmod;
    delete Jmod;
    delete JxB;
    delete st;
}

//-----------------------------------------------------------------------
bool PCOProcessor::setTask(ATQPTask * _task)
{
    bool finish = ATQPProcessor::setTask(_task);
    if (!finish)
        task = (PCOTask *)_task;

    return finish;
}

//-----------------------------------------------------------------------
double CagpWiegelmann::step(int _iterN, int _depth)
{
    iterN = _iterN;
    depth = _depth;
    counter++;

    out_field o(nProc, iterN, depth);
    o.out("vB", vB);

    mode = 0;
    main_proc->proceed(processors, supervisor, priority);

    o.out("B2", B2);
    o.out("invB", s);
    o.out("rotB", rotB);
    o.out("Wa", Wa);
    o.out("divB", divB);
    o.out("Wb", Wb);
    o.out("WaxB", WaxB);
    o.out("WbxB", WbxB);

    double s = 0;
    for (int kz = 0; kz < N[2]; kz++)
        s += L[kz];

    mode = 1;
    main_proc->proceed(processors, supervisor, priority);

    o.out("vF", (CagmVectorField *)vF);

    return s;
}

//-----------------------------------------------------------------------
bool PCOProcessor::proceed()
{
    if (host->mode == 0)
    {
        host->B2->abs2_plane(host->vB, task->z);
        host->B2sum[task->z] = host->B2->sum_plane(task->z);

        host->s->inv_plane(host->B2, task->z); // s1 == 1 / B^2

        // Wa
        host->rotB->rot_plane(host->vB, task->z, WiegelmannDerivStencil); // rotB ~ J
        host->Wa->cross_plane(host->rotB, host->vB, task->z); // rotB x B
        // Wb
        host->divB->div_plane(host->vB, task->z, WiegelmannDerivStencil);
        host->Wb->mult_plane(host->divB, host->vB, task->z); // divB * B
        if (WiegelmannGetMetricsTheta)
        {
            host->Bmod->sqrt_plane(host->B2, task->z);
            host->Binv->inv_plane(host->Bmod, task->z);
            host->Jmod->abs_plane(host->rotB, task->z);
            host->Jinv->inv_plane(host->Jmod, task->z);
            host->JxB->abs_plane(host->Wa, task->z);

            host->st->mult_plane(host->JxB, host->Binv, task->z);
            host->st->mult_plane(host->st, host->Jinv, task->z);
            host->sS[task->z] = host->st->sum_plane(task->z, host->sW);
            host->sSW[task->z] = host->sW->sum_plane(task->z);

            host->st->mult_plane(host->JxB, host->Binv, task->z);
            host->sSJ[task->z] = host->st->sum_plane(task->z, host->sW);
            host->sJ[task->z] = host->Jmod->sum_plane(task->z, host->sW);

            host->st->mult_plane(host->JxB, host->Jinv, task->z);
            host->sSB[task->z] = host->st->sum_plane(task->z, host->sW);
            host->sB[task->z] = host->Bmod->sum_plane(task->z, host->sW);
        }

        host->Wa->mult_plane(host->s, host->Wa, task->z); // rotB x B / B^2
        host->Wb->mult_plane(host->s, host->Wb, task->z); // divB * B / B^2

        // Wa^2 + Wb^2
        host->Wa2Wb2->abs2_plane(host->Wa, task->z);
        host->s->abs2_plane(host->Wb, task->z);
        if (WiegelmannGetMetricsTheta)
        {
            host->st->mult_plane(host->B2, host->Wa2Wb2, task->z);
            host->LF[task->z] = host->st->sum_plane(task->z, host->sW);
            host->st->mult_plane(host->B2, host->s, task->z);
            host->LD[task->z] = host->st->sum_plane(task->z, host->sW);
        }

        if (WiegelmannWeightDivfree != 1.0)
            host->s->mult_plane(WiegelmannWeightDivfree, host->s, task->z);
        host->Wa2Wb2->add_plane(host->Wa2Wb2, host->s, task->z);

        // functional
        host->s->mult_plane(host->B2, host->Wa2Wb2, task->z);

        if (task->z == 0)
            host->L[task->z] = 0;
        else
            host->L[task->z] = host->s->sum_plane(task->z, host->sW);

        host->WaxB->cross_plane(host->Wa, host->vB, task->z);
        host->WbxB->dot_plane(host->Wb, host->vB, task->z);
    }
    else
    {
        host->vF->mult_plane(host->Wa2Wb2, host->vB, task->z);
        host->v->rot_plane(host->WaxB, task->z, WiegelmannDerivStencil);
        host->vF->add_plane(host->vF, host->v, task->z);

        host->v->cross_plane(host->Wa, host->rotB, task->z);
        host->vF->sub_plane(host->vF, host->v, task->z);

        host->v->grad_plane(host->WbxB, task->z, WiegelmannDerivStencil);
        if (WiegelmannWeightDivfree != 1.0)
            host->v->mult_plane(WiegelmannWeightDivfree, host->v, task->z);
        host->vF->add_plane(host->vF, host->v, task->z); // rot(Wa x B) - Wa x rotB + grad(Wb . B) + (Wa^2 + Wb^2)*B

        host->v->mult_plane(host->divB, host->Wb, task->z);
        if (WiegelmannWeightDivfree != 1.0)
            host->v->mult_plane(WiegelmannWeightDivfree, host->v, task->z);
        host->vF->sub_plane(host->vF, host->v, task->z); //  rot(Wa x B) - Wa x rotB + grad(Wb . B) + (Wa^2 + Wb^2)*B - Wb*divB

        host->vF->mult_plane(host->sW, host->vF, task->z);

        // gradW
        host->v->cross_plane(host->WaxB, host->vgradW, task->z);
        host->vF->add_plane(host->vF, host->v, task->z);

        host->v->mult_plane(host->WbxB, host->vgradW, task->z);
        if (WiegelmannWeightDivfree != 1.0)
            host->v->mult_plane(WiegelmannWeightDivfree, host->v, task->z);
        host->vF->add_plane(host->vF, host->v, task->z); // rot(Wa x B) x gradW + (Wb . B)*gradW
        host->F2max[task->z] = host->vF->max2_plane(task->z);
    }

    return true;
}
