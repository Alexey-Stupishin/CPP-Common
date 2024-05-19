#include "TaskQueueProcessor.h"

//-----------------------------------------------------------------------------
ATQPTask::ATQPTask(size_t task_id)
    : id(task_id) 
{
    console_debug("task " << id << " created")
}

//-----------------------------------------------------------------------------
ATQPTask::~ATQPTask() 
{
    console_debug("task " << id << " deleted")
}

//-----------------------------------------------------------------------------
size_t ATQPTask::getID() 
{ 
    return id; 
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
ATQPTaskFactory::ATQPTaskFactory()
    : counter(0)
{}

//-----------------------------------------------------------------------------
ATQPTaskFactory::~ATQPTaskFactory()
{}

////-----------------------------------------------------------------------------
//void ATQPTaskFactory::initialize(ATQPTask *task)
//{ 
//    task->id = 11 + counter++; 
//} // ToDo: can be some uuid generation
