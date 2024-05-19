#ifdef USE_UNIT_TEST

/*
subclassing

class QPTask : public ATQPTask
    void setData(int _data);
    int getData();

class QPTaskFactory : public ATQPTaskFactory
    virtual ATQPTask *create()

class QPSupervisor : public ATQPSupervisor
    QPSupervisor

class QPProcessor : public ATQPProcessor
    virtual bool setTask(ATQPTask * _task)
    virtual bool proceed()
*/

/*
Sample console output for call TaskQueueProcessor::unitTest(2, 3):

sync objects cleared
supervisor objects cleared
started
sync objects cleared
supervisor objects cleared
[supervisor]    running...
[processor 0]   running...
[processor 1]   running...
[supervisor]    query task from processor 0, task address 000002BE961A5E00
  [child processor 0]   get task 0 (task address 000002BE961A5E00), task data 0
[processor 0]   proceeds task by address 000002BE961A5E00
  [child processor 1]   get task 1 (task address 000002BE961A63A0), task data 2
[processor 1]   proceeds task by address 000002BE961A63A0
  [child processor 1]   task 1 will be proceeded in 4 seconds
[supervisor]    query task from processor 1, task address 000002BE961A63A0
  [child processor 0]   task 0 will be proceeded in 4 seconds
[supervisor]    query task from processor 1, task address 000002BE961A5800
  [child processor 1]   get task 2 (task address 000002BE961A5800), task data 4
[processor 1]   proceeds task by address 000002BE961A5800
  [child processor 1]   task 2 will be proceeded in 4 seconds
[supervisor]    query task from processor 0, task address 0000000000000000
[processor 0]   remain(s) 1 active processor(s)
[supervisor]    query task from processor 1, task address 0000000000000000
[processor 1]   remain(s) 0 active processor(s)
supervisor done
end of supervisor
terminated 2 processors
finished
sync objects deleted
total finish
*/

#include "TaskQueueProcessor.h"

#include <random>
#include <chrono>


//---------------------------------------------------------------------------------------
class QPTask : public ATQPTask
{
protected:
    int data;

public:
    QPTask(int _id) : ATQPTask(_id) {}
    void setData(int _data) { data = _data; }
    int getData() { return data; }
};
    
//---------------------------------------------------------------------------------------
class QPTaskFactory : public ATQPTaskFactory
{
public:
    QPTaskFactory() : ATQPTaskFactory() {}
    virtual ~QPTaskFactory() {}

    virtual ATQPTask *create()
    {
        return new QPTask(counter++);
        //initialize(task);

        //return task;
    }
};

//---------------------------------------------------------------------------------------
class QPSupervisor : public ATQPSupervisor
{
public:
    QPSupervisor(int _n_task, ATQPTaskFactory* _factory, ATQPSynchonizer *_sync) 
        : ATQPSupervisor(_n_task, _factory, _sync)
    {
        for (int i = 0; i < _n_task; i++)
            ((QPTask *)tasks[i])->setData(2*i);
    }

    virtual ~QPSupervisor() 
    { 
    }
};

//---------------------------------------------------------------------------------------
class QPProcessor : public ATQPProcessor
{
protected:
    std::mt19937 *p_generator;
    QPTask *this_task;

public:
    QPProcessor(int _id, ATQPSynchonizer *_sync) : ATQPProcessor(_id, _sync)
    {
        p_generator = new std::mt19937((unsigned int)std::chrono::system_clock::now().time_since_epoch().count());
    }
    virtual ~QPProcessor() { delete p_generator; }

    virtual bool setTask(ATQPTask * _task)
    {
        bool finish = ATQPProcessor::setTask(_task);
        if (!finish)
        {   // cast void *_task to the real object:
            this_task = (QPTask *)_task;
            console_debug("  [child processor " << id << "]\tget task " << this_task->getID() << " (task address " << _task << "), task data " << this_task->getData());
        }
        return finish;
    }

    virtual bool proceed()
    {
        int delay = this_task->getData() /2 + (*p_generator)() % 5;
        console_debug("  [child processor " << id << "]\ttask " << this_task->getID() << " will be proceeded in " << delay << " seconds");
        std::this_thread::sleep_for(std::chrono::seconds(delay));
        return true;
    }
};

//---------------------------------------------------------------------------------------
void TaskQueueProcessor::unitTest2(int num_proc, int num_tasks, int num_loops)
{
    console_start();

    TaskQueueProcessor proc(num_proc);

    std::vector<ATQPProcessor *> processors;
    for (int i = 0; i < num_proc; i++)
        processors.push_back(new QPProcessor(i, proc.get_sync()));

    QPTaskFactory factory;
    QPSupervisor supervisor(num_tasks, &factory, proc.get_sync());
    for (int k = 0; k < num_loops; k++)
    {
        console_debug("started");
        proc.proceed(processors, &supervisor);
        console_debug("finished");
    }

    for (int i = 0; i < num_proc; i++)
        delete processors[i];
}

//---------------------------------------------------------------------------------------
void TaskQueueProcessor::unitTest(int num_proc, int num_tasks)
{
    console_start();

    TaskQueueProcessor *proc = new TaskQueueProcessor(num_proc);

    std::vector<ATQPProcessor *> processors;
    for (int i = 0; i < num_proc; i++)
        processors.push_back(new QPProcessor(i, proc->get_sync()));

    QPTaskFactory factory;
    QPSupervisor *supervisor = new QPSupervisor(num_tasks, &factory, proc->get_sync());

    console_debug("started");
    proc->proceed(processors, supervisor);
    console_debug("finished");

    for (int i = 0; i < num_proc; i++)
        delete processors[i];

    delete supervisor;
    delete proc;

    console_debug("total finish");
}

#endif // USE_UNIT_TEST
