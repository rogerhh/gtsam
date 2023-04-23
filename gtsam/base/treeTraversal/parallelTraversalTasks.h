/* ----------------------------------------------------------------------------

* GTSAM Copyright 2010, Georgia Tech Research Corporation,
* Atlanta, Georgia 30332-0415
* All Rights Reserved
* Authors: Frank Dellaert, et al. (see THANKS for the full author list)

* See LICENSE for the license information

* -------------------------------------------------------------------------- */

/**
* @file    parallelTraversalTasks.h
* @author  Richard Roberts
* @date    April 9, 2013
*/
#pragma once

#include <gtsam/global_includes.h>

#include <memory>

#ifdef GTSAM_USE_TBB
#include <tbb/task_group.h>         // tbb::task_group
#include <tbb/scalable_allocator.h> // tbb::scalable_allocator
// #include <tbb/tbb_thread.h>
#include <tbb/task_arena.h>
#include <iostream>

namespace gtsam {

  /** Internal functions used for traversing trees */
  namespace treeTraversal {

    namespace internal {

      /* ************************************************************************* */
      template<typename NODE, typename DATA, typename VISITOR_PRE, typename VISITOR_POST>
      class PreOrderTask
      {
      public:
        const std::shared_ptr<NODE>& treeNode;
        std::shared_ptr<DATA> myData;
        VISITOR_PRE& visitorPre;
        VISITOR_POST& visitorPost;
        int problemSizeThreshold;
        tbb::task_group& tg;
        bool makeNewTasks;

        // Keep track of order phase across multiple calls to the same functor
        mutable bool isPostOrderPhase;

        PreOrderTask(const std::shared_ptr<NODE>& treeNode, const std::shared_ptr<DATA>& myData,
                     VISITOR_PRE& visitorPre, VISITOR_POST& visitorPost, int problemSizeThreshold,
                     tbb::task_group& tg, bool makeNewTasks = true)
            : treeNode(treeNode),
              myData(myData),
              visitorPre(visitorPre),
              visitorPost(visitorPost),
              problemSizeThreshold(problemSizeThreshold),
              tg(tg),
              makeNewTasks(makeNewTasks),
              isPostOrderPhase(false) {}

        void operator()() const
        {
          if(isPostOrderPhase)
          {
            // Run the post-order visitor since this task was recycled to run the post-order visitor
            (void) visitorPost(treeNode, *myData);
          }
          else
          {
            if(makeNewTasks)
            {
              if(!treeNode->children.empty())
              {
                bool overThreshold = (treeNode->problemSize() >= problemSizeThreshold);

                // If we have child tasks, start subtasks and wait for them to complete
                tbb::task_group ctg;
                fflush(stdout);
                printf("ARENA: %d; thread: %p, Parent task group INNER: %p, Child task group START: %p, size: %zu\n",
                    tbb::this_task_arena::current_thread_index(),
                    std::this_thread::get_id(),
                    &tg,
                    &ctg,
                    treeNode->children.size());
                // std::cout 
                //     << "ARENA: " << tbb::this_task_arena::current_thread_index()
                //     << ": "
                //     << "thread " << std::this_thread::get_id()
                //     << ", Child task group START " << &ctg
                //     << ", size " << (int)treeNode->children.size()
                //     << "\n"
                //     << std::endl;
                fflush(stdout);

                for(const std::shared_ptr<NODE>& child: treeNode->children)
                {
                  // Process child in a subtask.  Important:  Run visitorPre before calling
                  // allocate_child so that if visitorPre throws an exception, we will not have
                  // allocated an extra child, this causes a TBB error.
                  std::shared_ptr<DATA> childData = std::allocate_shared<DATA>(
                      tbb::scalable_allocator<DATA>(), visitorPre(child, *myData));
                  ctg.run(PreOrderTask(child, childData, visitorPre, visitorPost,
                      problemSizeThreshold, ctg, overThreshold));
                }
                ctg.wait();

                fflush(stdout);
                printf("ARENA: %d; thread: %p, Parent task group INNER: %p, Child task group END: %p, size: %zu\n",
                    tbb::this_task_arena::current_thread_index(),
                    std::this_thread::get_id(),
                    &tg,
                    &ctg,
                    treeNode->children.size());
                // std::cout 
                //     << "ARENA: " << tbb::this_task_arena::current_thread_index()
                //     << ": "
                //     << "thread " << std::this_thread::get_id()
                //     << ", Child task group END " << &ctg
                //     << ", size " << (int)treeNode->children.size()
                //     << "\n"
                //     << std::endl;
                fflush(stdout);

                // Allocate post-order task as a continuation
                isPostOrderPhase = true;
                tg.run(*this);
              }
              else
              {

                // Run the post-order visitor in this task if we have no children
                fflush(stdout);
                printf("ARENA: %d; thread: %p, Parent task group LEAF: %p\n",
                    tbb::this_task_arena::current_thread_index(),
                    std::this_thread::get_id(),
                    &tg);
                // std::cout 
                //     << "ARENA: " << tbb::this_task_arena::current_thread_index()
                //     << ": "
                //     << "thread " << std::this_thread::get_id()
                //     << ", Parent task group LEAF " << &tg
                //     << "\n"
                //     << std::endl;
                fflush(stdout);
                (void) visitorPost(treeNode, *myData);
              }
            }
            else
            {
              // Process this node and its children in this task
              processNodeRecursively(treeNode, *myData);
            }
          }
        }

        void processNodeRecursively(const std::shared_ptr<NODE>& node, DATA& myData) const
        {
          for(const std::shared_ptr<NODE>& child: node->children)
          {
            DATA childData = visitorPre(child, myData);
            processNodeRecursively(child, childData);
          }

          // Run the post-order visitor
          (void) visitorPost(node, myData);
        }
      };

      /* ************************************************************************* */
      template<typename ROOTS, typename NODE, typename DATA, typename VISITOR_PRE, typename VISITOR_POST>
      class RootTask
      {
      public:
        const ROOTS& roots;
        DATA& myData;
        VISITOR_PRE& visitorPre;
        VISITOR_POST& visitorPost;
        int problemSizeThreshold;
        tbb::task_group& tg;
        RootTask(const ROOTS& roots, DATA& myData, VISITOR_PRE& visitorPre, VISITOR_POST& visitorPost,
          int problemSizeThreshold, tbb::task_group& tg) :
          roots(roots), myData(myData), visitorPre(visitorPre), visitorPost(visitorPost),
          problemSizeThreshold(problemSizeThreshold), tg(tg) {}

        void operator()() const
        {
          typedef PreOrderTask<NODE, DATA, VISITOR_PRE, VISITOR_POST> PreOrderTask;
          // Create data and tasks for our children
          for(const std::shared_ptr<NODE>& root: roots)
          {
            std::shared_ptr<DATA> rootData = std::allocate_shared<DATA>(tbb::scalable_allocator<DATA>(), visitorPre(root, myData));
            tg.run(PreOrderTask(root, rootData, visitorPre, visitorPost, problemSizeThreshold, tg));
          }
        }
      };

      template<typename NODE, typename ROOTS, typename DATA, typename VISITOR_PRE, typename VISITOR_POST>
      void CreateRootTask(const ROOTS& roots, DATA& rootData, VISITOR_PRE& visitorPre, VISITOR_POST& visitorPost, int problemSizeThreshold)
      {
          typedef RootTask<ROOTS, NODE, DATA, VISITOR_PRE, VISITOR_POST> RootTask;
          tbb::task_group tg;
          printf("Root Task created with %zu roots\n", roots.size());
          tg.run_and_wait(RootTask(roots, rootData, visitorPre, visitorPost, problemSizeThreshold, tg));
          printf("Root Task finished with %zu roots\n\n\n", roots.size());
      }

    }

  }

}

#endif
