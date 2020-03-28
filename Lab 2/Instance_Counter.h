#pragma once

#ifndef _Instance_Counter_h
#define _Instance_Counter_h

class Instance_Counter 
{
    public:
        static int counter;
        
        Instance_Counter(void)
        {
            counter++;
        }
};

#endif /* _Instance_Counter_h */