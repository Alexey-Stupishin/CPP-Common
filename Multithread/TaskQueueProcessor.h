#pragma once

#include <vector>
#include <cstddef>

#include <thread>
#include <queue>
#include <condition_variable>
#include <mutex>

#include "console_debug.h"

#ifdef _WINDOWS
#include <windows.h>
#endif

enum w_priority /* : int */ { lowest = -2, low = -1, normal = 0, high = 1, highest = 2 };

class ATQPTask;

//---------------------------------------------------------------------------------------
class ATQPSynchonizer
{
protected:
    std::size_t num_proc;

public:
    std::size_t             process_counter; // processor, supervisor, TQP.proceed

    std::mutex              mutex_query; // processor, supervisor
    std::condition_variable check_query; // processor, supervisor
    std::queue<int>         queue_query; // processor, supervisor, TQP.proceed

    std::vector<std::mutex *>              mutexes_task; // processor, supervisor, TQP.proceed
    std::vector<std::condition_variable *> checks_task; // processor, supervisor, TQP.proceed
    std::vector<std::queue<ATQPTask *> *>  queues_task; // processor, supervisor, TQP.proceed

    ATQPSynchonizer(std::size_t num_proc);
    virtual ~ATQPSynchonizer();
    void reset();
};

//---------------------------------------------------------------------------------------
class ATQPTask
{
friend class ATQPTaskFactory;

protected:
    size_t id;

public:
    ATQPTask(size_t task_id = 0);
    virtual ~ATQPTask();

    size_t getID();
};

//---------------------------------------------------------------------------------------
class ATQPTaskFactory
{
protected:
    int counter;

public:
    ATQPTaskFactory();
    virtual ~ATQPTaskFactory();

    virtual ATQPTask *create() = 0;
    //void initialize(ATQPTask *task);
};

//---------------------------------------------------------------------------------------
class ATQPSupervisor
{
protected:
    ATQPTask **tasks;
    int n_task;
    int current_task;
    ATQPTaskFactory *factory;

public:
    ATQPSynchonizer *sync;

public:
    ATQPSupervisor(int _n_task, ATQPTaskFactory *_factory, ATQPSynchonizer *_sync);
    virtual ~ATQPSupervisor();

    bool reset();
    virtual bool getTask(ATQPTask*& task);
};

//---------------------------------------------------------------------------------------
class ATQPProcessor
{
protected:
    ATQPTask *task;
    int id;
    bool init;

public:
    ATQPSynchonizer *sync;

public:
    ATQPProcessor(int _id, ATQPSynchonizer *_sync);
    virtual ~ATQPProcessor();

    bool initialized();
    bool reset();
    int getID();
    virtual bool setTask(ATQPTask *_task);

    virtual ATQPTask* getTask();
    virtual bool proceed() = 0;
};

//---------------------------------------------------------------------------------------
class TaskQueueProcessor
{
protected:
    ATQPSynchonizer *sync;
    int num_proc;

public:
    TaskQueueProcessor(int nThreadsInitial = 0);
    virtual ~TaskQueueProcessor();

    unsigned long proceed(std::vector<ATQPProcessor *>&, ATQPSupervisor *, w_priority priority);
    ATQPSynchonizer *get_sync();
    int get_num_proc();

    static int getProcInfo(int nThreadsInitial = 0);

#ifdef USE_UNIT_TEST
    static void unitTest(int num_proc = 5, int num_tasks = 10);
    static void unitTest2(int num_proc = 2, int num_tasks = 3, int num_loops = 5);
#endif
};
