#include "TaskQueueProcessor.h"

//#include <thread>
//#include "console_debug.h"

//-----------------------------------------------------------------------------
ATQPSynchonizer::ATQPSynchonizer(std::size_t _num_proc)
{
    num_proc = _num_proc;
    for (std::size_t i = 0; i < num_proc; i++)
    {
        mutexes_task.push_back(new std::mutex);
        checks_task.push_back(new std::condition_variable);
        queues_task.push_back(new std::queue<ATQPTask *>);
    }

    console_debug("sync object created, num_proc = " << num_proc)

    reset();
}

//-----------------------------------------------------------------------------
ATQPSynchonizer::~ATQPSynchonizer()
{
    for (int i = 0; i < num_proc; i++)
    {
        delete mutexes_task[i];
        delete checks_task[i];
        delete queues_task[i];
    }

    console_debug("sync object deleted")
}

//-----------------------------------------------------------------------------
void ATQPSynchonizer::reset()
{
    queue_query.empty();
    process_counter = num_proc;

    console_debug("sync object reset")
}
