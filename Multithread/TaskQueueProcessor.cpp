#include "TaskQueueProcessor.h"

//-----------------------------------------------------------------------------
void processorFunc(ATQPProcessor *processor)
{
    int id = processor->getID();
    console_debug("[processor " << id << "]\trunning...")

    bool finish = false;
    while (!finish)
    {
        if (processor->initialized())
        {
            console_debug("[processor " << id << "]\tproceeds task by address " << processor->getTask());
            processor->proceed();
        }

        { // supervisor request
            std::unique_lock<std::mutex> locker(processor->sync->mutex_query);
            processor->sync->queue_query.push(id);
            processor->sync->check_query.notify_one();
        }

        { // supervisor respond
            std::unique_lock<std::mutex> locker(*(processor->sync->mutexes_task[id]));
            processor->sync->checks_task[id]->wait(locker, [&]() {return !processor->sync->queues_task[id]->empty();});
            if (!processor->sync->queues_task[id]->empty())
            {
                finish = processor->setTask(processor->sync->queues_task[id]->front());
                processor->sync->queues_task[id]->pop();
            }
        }
    }

    {
        std::unique_lock<std::mutex> locker(processor->sync->mutex_query);
        processor->sync->process_counter--;
        console_debug("[processor " << id << "]\tremain(s) " << processor->sync->process_counter << " active processor(s)");
        processor->reset();
        processor->sync->check_query.notify_one();
    }

    console_debug("end of processor thread");
}

//-----------------------------------------------------------------------------
void supervisorFunc(ATQPSupervisor* supervisor)
{
    console_debug("[supervisor]\trunning...");

    bool done = false;
    while (!done)
    {
        std::unique_lock<std::mutex> locker(supervisor->sync->mutex_query);
        supervisor->sync->check_query.wait(locker, [&]() { return !supervisor->sync->queue_query.empty() || supervisor->sync->process_counter == 0; });
        if (supervisor->sync->process_counter == 0)
        {
            done = true;
            console_debug("supervisor done");
        }
        while (!supervisor->sync->queue_query.empty() && !done)
        {
            int id = supervisor->sync->queue_query.front();
            ATQPTask* task;
            {
                std::unique_lock<std::mutex> locker_m(*(supervisor->sync->mutexes_task[id]));
                supervisor->getTask(task);
                supervisor->sync->queues_task[id]->push(task);
                supervisor->sync->checks_task[id]->notify_one();
            }
            supervisor->sync->queue_query.pop();
            if (task)
                console_debug("[supervisor]\tquery task " << task->getID() << " from processor " << id << ", task address " << task)
            else
                console_debug("[supervisor]\tprocessor " << id << ", no task")
        }
    }

    console_debug("end of supervisor thread");
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
TaskQueueProcessor::TaskQueueProcessor(int nThreadsInitial)
{
    num_proc = TaskQueueProcessor::getProcInfo(nThreadsInitial);
    sync = new ATQPSynchonizer(num_proc);
}

//-----------------------------------------------------------------------------
TaskQueueProcessor::~TaskQueueProcessor()
{
    delete sync;

    console_debug("main processor object deleted")

}

//-----------------------------------------------------------------------------
unsigned long TaskQueueProcessor::proceed(std::vector<ATQPProcessor *>& processors, ATQPSupervisor *supervisor, w_priority priority)
{

    std::size_t _num_proc = processors.size();

    sync->reset();
    supervisor->reset();

    std::thread supervisorThread(supervisorFunc, supervisor);
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < _num_proc; i++)
        threads.push_back(std::thread(processorFunc, processors[i]));

    for (auto &t : threads)
    {
#ifdef _WINDOWS
        int sys_p;
        switch (priority)
        {
            case lowest:
                sys_p = THREAD_PRIORITY_LOWEST;
                break;
            case low:
                sys_p = THREAD_PRIORITY_BELOW_NORMAL;
                break;
            case normal:
                sys_p = THREAD_PRIORITY_NORMAL;
                break;
            case high:
                sys_p = THREAD_PRIORITY_ABOVE_NORMAL;
                break;
            case highest:
                sys_p = THREAD_PRIORITY_HIGHEST;
                break;
            default:
                sys_p = THREAD_PRIORITY_BELOW_NORMAL;
        }
        HANDLE h = (HANDLE)(t.native_handle());
        SetThreadPriority(h, sys_p);
#endif
        t.join();
    }
    supervisorThread.join();

    console_debug("terminated " << _num_proc << " processors")

    return 0;
}

//-----------------------------------------------------------------------------
ATQPSynchonizer * TaskQueueProcessor::get_sync()
{
    return sync;
}

//-----------------------------------------------------------------------------
int TaskQueueProcessor::get_num_proc()
{ 
    return num_proc;  
}

//-----------------------------------------------------------------------------
int TaskQueueProcessor::getProcInfo(int nThreadsInitial)
{
    int n = (int)std::thread::hardware_concurrency();

    int outn = 1;
    if (nThreadsInitial > 100000) // debug!
        outn = nThreadsInitial - 100000;
    else if (nThreadsInitial < 0)
        outn = n - nThreadsInitial;
    else if (nThreadsInitial == 0)
        outn = n;
    else
    {
        if (nThreadsInitial > n)
            outn = n;
        else
            outn = nThreadsInitial;
    }

    if (outn <= 0)
        outn = 1;

    return outn;
}
