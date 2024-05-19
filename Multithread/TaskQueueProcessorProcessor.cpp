#include "TaskQueueProcessor.h"

//-----------------------------------------------------------------------------
ATQPProcessor::ATQPProcessor(int _id, ATQPSynchonizer *_sync) 
        : id(_id)
        , init(false)
        , sync(_sync)
{}

//-----------------------------------------------------------------------------
ATQPProcessor::~ATQPProcessor() 
{}

//-----------------------------------------------------------------------------
bool ATQPProcessor::initialized() 
{ 
    return init; 
}

//-----------------------------------------------------------------------------
bool ATQPProcessor::reset()
{ 
    init = false; 
    return init; 
}

//-----------------------------------------------------------------------------
int ATQPProcessor::getID()
{ 
    return id; 
}

//-----------------------------------------------------------------------------
bool ATQPProcessor::setTask(ATQPTask *_task)
{
    if (!_task)
        return true;
    task = _task;
    init = true;
    return false;
}

//-----------------------------------------------------------------------------
ATQPTask* ATQPProcessor::getTask() 
{ 
    return task; 
}
