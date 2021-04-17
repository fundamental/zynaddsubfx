/*
  ZynAddSubFX - a software synthesizer

  MultiPseudoStack.h - Multiple-Writer Lock Free Datastructure
  Copyright (C) 2016 Mark McCurry

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/
#pragma once
#include <atomic>

namespace zyn {

//XXX rename this thing
typedef struct QueueListItem qli_t;
//! Queue item for LockFreeQueue
struct QueueListItem
{
    QueueListItem(void);
    char     *memory;
    uint32_t  size;
};


//! Many reader many writer queue of fixed size memory chunks
class LockFreeQueue
{
    qli_t *const           data;
    const int              elms;
    std::atomic<uint32_t> *tag;
    std::atomic<int32_t>   next_r;
    std::atomic<int32_t>   next_w;
    std::atomic<int32_t>   avail;
public:
    LockFreeQueue(qli_t *data_, int n);
    ~LockFreeQueue(void);
    qli_t *read(void);
    void write(qli_t *Q);
    int8_t available();
};


/**
 * Many reader many writer capiable queue of fixed size memory chunks
 *
 * - lock free
 * - allocation free (post initialization)
 */
class MultiQueue
{
    qli_t *pool;          //!< actual storage
    LockFreeQueue m_free; //!< queue of unused items
    LockFreeQueue m_msgs; //!< actual queue

public:
    MultiQueue(void);
    ~MultiQueue(void);
    //! return a fresh item, e.g. to write to
    qli_t *alloc(void)   { return m_free.read();   }
    //! free an item
    void free(qli_t  *q) {        m_free.write(q); }
    //! write item into the queue
    void write(qli_t *q) {        m_msgs.write(q); }
    //! take an item from the queue, must be freed
    qli_t *read(void)    { return m_msgs.read();   }
    //! return if there's at least one item in the queue
    bool canRead()       { return m_msgs.available(); }
};

}
