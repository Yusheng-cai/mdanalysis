#pragma once
#include "OpenMP.h"
#include <vector>

namespace OpenMP{
template<typename T>
class OpenMP_buffer{
    public:
        // iterator & const iterator of the vector

        // Explanation for adding typename 
        // The compiler seems to be saying that's it's syntactically unclear 
        // (at that point in the compilation anyway) whether or not vector<T>::iterator is a typename
        using iterator        = typename std::vector<T>::iterator;
        using const_iterator  = typename std::vector<T>::const_iterator;

        // The number of usable threads is total number of threads
        OpenMP_buffer(){buffer_.resize(OpenMP::get_max_threads());}
        OpenMP_buffer(T& master_obj):master_object_ptr_(&master_obj){buffer_.resize(OpenMP::get_max_threads());}

        // set the mater object
        void set_master_object(T& master_obj){master_object_ptr_ = &master_obj;}

        // Access the specific buffer by thread id
        T& access_buffer_by_id(){
            int id = OpenMP::get_thread_num();
            if (id == 0 && master_object_ptr_!= nullptr){
                // we are returning the reference
                // pointer and reference are different 
                // pointer has its own memory address but reference shares the mem address with the original obj
                return *master_object_ptr_;
            }
            else{
                return buffer_[id];
            }
        }

        // clear the master_object_ptr
        void clearMasterObject(){master_object_ptr_=nullptr;}

        // clear all buffers
        void clearBuffer()
        {
            if (buffer_.size() != 0)
            {
                for (int i=0;i<buffer_.size();i++)
                {
                    buffer_[i].clear();
                }
            }
        }

        T* getMasterObject(){return master_object_ptr_;}

        // start from the second buffer since the first one is always the "master object"
        iterator beginworker() {return buffer_.begin()++;}
        iterator endworker()   {return buffer_.end();}

    private:
        std::vector<T> buffer_; 
        T*  master_object_ptr_ = nullptr;
};
}