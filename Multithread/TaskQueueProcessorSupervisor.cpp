#include "TaskQueueProcessor.h"

//-----------------------------------------------------------------------------
ATQPSupervisor::ATQPSupervisor(int _n_task, ATQPTaskFactory *_factory, ATQPSynchonizer *_sync) 
    : n_task(_n_task)
    , tasks(nullptr)
    , factory(_factory)
    , sync(_sync)
{ 
    if (n_task > 0)
    {
        tasks = new ATQPTask *[n_task];
        for (int i = 0; i < n_task; i++)
            tasks[i] = factory->create();
    }

    reset();
}

//-----------------------------------------------------------------------------
ATQPSupervisor::~ATQPSupervisor()
{
    for (int i = 0; i < n_task; i++)
        delete tasks[i];

    delete[] tasks;

    console_debug("supervisor objects deleted")
}

//-----------------------------------------------------------------------------
bool ATQPSupervisor::ATQPSupervisor::reset()
{
    current_task = 0;

    console_debug("supervisor objects reset")

    return true;
}

//-----------------------------------------------------------------------------
bool ATQPSupervisor::getTask(ATQPTask*& task)
{
    task = nullptr;
    if (current_task >= n_task)
        return false;
    task = tasks[current_task];
    current_task++;
    return true;
}
