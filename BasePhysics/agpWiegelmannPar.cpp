#include "stdDefinitions.h"
#include <math.h>

#include "agpWiegelmannPar.h"
#include "agmScalarField.h"
#include "agmVectorField.h"
#include "mfoGlobals.h"

#include "DebugWrite.h"


//-----------------------------------------------------------------------
CagpWiegelmann::CagpWiegelmann(int *_N, int _n_threads, int _stencil
    , CagmVectorField *_sourceB, CagmScalarField *_sourceW
    , double *_vcos
    , CagmVectorFieldOps *_outF
    , CagmVectorField *_baseField, CagmVectorField *_baseWeight
    , CagmScalarField *_absField, CagmScalarField *_absWeight
    , CagmScalarField *_losField, CagmScalarField *_losWeight
    , CagmScalarField *_bottomWeight
    , int _depth
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
    , depth(_depth)
    , priority(_priority)
{
    vB->dimensions(N);

    vgradW = new CagmVectorField(sW);

    DebugWriteData(sW, "debug_W", depth);

    vgradW->grad(sW);

    DebugWriteData(vgradW, "debug_gradW", depth);

    int Nb[3];
    Nb[0] = N[0];
    Nb[1] = N[1];
    Nb[2] = 1;
    vBottom = new CagmVectorField(Nb);

    B2 = new CagmScalarField(vB);
    rotB = new CagmVectorField(vB);
    divB = new CagmScalarField(vB);
    Wa = new CagmVectorField(vB);
    Wb = new CagmVectorField(vB);
    Wa2Wb2 = new CagmScalarField(vB);
    WaxB = new CagmVectorField(vB);
    WbxB = new CagmScalarField(vB);
    v = new CagmVectorField(vB);
    s = new CagmScalarField(vB);

    if (WiegelmannGetMetricsTheta)
    {
        st  =  new CagmScalarField(vB);
        Binv = new CagmScalarField(vB);
        Jinv = new CagmScalarField(vB);
        Bmod = new CagmScalarField(vB);
        Jmod = new CagmScalarField(vB);
        JxB  = new CagmScalarField(vB);
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
double CagpWiegelmann::step(int _iterN)
{
    iterN = _iterN;
    counter++;

    //DebugWriteData(vB, "vB", depth, iterN);

    mode = 0;
    main_proc->proceed(processors, supervisor, priority);

    if (depth == 1 && iterN == 0)
    { 
        DebugWriteData(B2, "B2", depth, iterN);
        //DebugWriteData(invB, "invB", depth, iterN);
        DebugWriteData(rotB, "rotB", depth, iterN);
        DebugWriteData(Wa, "Wa", depth, iterN);
        DebugWriteData(divB, "divB", depth, iterN);
        DebugWriteData(Wb, "Wb", depth, iterN);
        DebugWriteData(WaxB, "WaxB", depth, iterN);
        DebugWriteData(WbxB, "WbxB", depth, iterN);
    }

    double _s = 0;
    for (int kz = 0; kz < N[2]; kz++)
        _s += L[kz];

    mode = 1;
    main_proc->proceed(processors, supervisor, priority);

    //DebugWriteData(vF, "vF", depth, iterN);

    return _s;
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
