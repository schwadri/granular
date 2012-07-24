/** \file time.cpp
 *  \author Adrian Schweizer
 *  \created  $Sa 18 Aug 06:32:47 pm CEST 2007 schwadri@SchwadriComp.local$
 *  \modified $So 30 Sep 04:44:49 pm CEST 2007 schwadri@SchwadriLaptop$
 *  \par FIXME
 *  - msvc++ does not support c99 and therefore uintxx_t is not defined (yepee! great)
 *    find alternative solution
 */

#include "time.hpp"
#include <limits>

//@{from boost::thread xtime.cpp (slightly changed)

#if defined(BOOST_HAS_FTIME)
#   define __STDC_CONSTANT_MACROS
#endif

#if defined(BOOST_HAS_FTIME)
#   include <windows.h>
#   include <boost/cstdint.hpp>
//for some stupid reason the windows header file define a macro called max. good job
#undef max
#elif defined(BOOST_HAS_GETTIMEOFDAY)
#   include <sys/time.h>
#elif defined(BOOST_HAS_MPTASKS)
#error "Support for anything prior MacOsX 1.0 is deprecated"
#endif

#include <cassert>
//@}


namespace core {

    const   time    time::PAST      =   time(0,0);
    const   time    time::FUTURE    =   time(std::numeric_limits<time::sec_type>::max(),std::numeric_limits<time::nsec_type>::max() );
    int const time::NSEC_PER_SEC   =   1000000000;
    int const time::NSEC_PER_MSEC  =   1000000;
    const   int time::NSEC_PER_USEC  =   1000;
    const   int time::USEC_PER_MSEC  =   1000;
    const   int time::MSEC_PER_SEC   =   1000;
    const   int time::USEC_PER_SEC   =   1000000;


    time    operator+(const time& lt, const time& rt)
    {
        time    temp(lt);
        temp    +=  rt;
        return temp;
    }

    time    operator-(const time& lt, const time& rt)
    {
        time    temp(lt);
        temp    -=  rt;
        return temp;
    }

    bool    operator==(const time& lt, const time& rt)
    {
        return  lt.nsec == rt.nsec && lt.sec == rt.sec;
    }

    bool    operator!=(const time& lt, const time& rt)
    {
        return  !(operator==(lt,rt));
    }

    bool    operator<(const time& lt, const time& rt)
    {
        return  lt.sec < rt.sec || (lt.sec == rt.sec && lt.nsec< rt.nsec);
    }

    bool    operator>(const time& lt, const time& rt)
    {
        return  lt.sec > rt.sec || (lt.sec == rt.sec && lt.nsec > rt.nsec);
    }

    bool    operator<=(const time& lt, const time& rt)
    {
        return  operator<(lt,rt) || operator==(lt,rt);
    }

    bool    operator>=(const time& lt, const time& rt)
    {
        return  operator>(lt,rt) || operator==(lt,rt);
    }

    void time::get(time& xtp)
    {
        
    //TODO: this is just one possible solution for the windows platform
    //      check whether another one would be faster/more precise
    #if defined(BOOST_HAS_FTIME)
		static bool s_initialized=false;
    /*  performance counter has better resolution...
        FILETIME ft;
    #if defined(BOOST_NO_GETSYSTEMTIMEASFILETIME)
        {
            SYSTEMTIME st;
            GetSystemtime(&st);
            SystemtimeToFiletime(&st,&ft);
        }
    #else
        GetSystemtimeAsFiletime(&ft);
    #endif
        static const boost::uint64_t TIMESPEC_TO_FILETIME_OFFSET = UINT64_C(116444736000000000);

        const boost::uint64_t ft64 = (static_cast<boost::uint64_t>(ft.dwHighDatetime) << 32)
            + ft.dwLowDatetime;

        xtp.sec = static_cast<xtime::xtime_sec_t>((ft64 - TIMESPEC_TO_FILETIME_OFFSET) / 10000000);

        xtp.nsec = static_cast<xtime::xtime_nsec_t>(((ft64 - TIMESPEC_TO_FILETIME_OFFSET) % 10000000) * 100);
    */

    static LARGE_INTEGER s_init_time;
    static LARGE_INTEGER s_to_nsec;
    static LARGE_INTEGER s_to_sec;
    if(!s_initialized)
    {
        LARGE_INTEGER c_freq;
        QueryPerformanceCounter(&s_init_time);
        QueryPerformanceFrequency(&c_freq);
        //prepare this value so that division through it directly leads to milliseconds
        s_to_nsec.QuadPart  = c_freq.QuadPart/NSEC_PER_SEC;
        s_to_sec.QuadPart   = c_freq.QuadPart;
        s_initialized = true;
    }
    LARGE_INTEGER t;
    QueryPerformanceCounter(&t);
    xtp.sec     = static_cast<nsec_type>(t.QuadPart/s_to_sec.QuadPart);
    xtp.nsec    = static_cast<nsec_type>(t.QuadPart/s_to_nsec.QuadPart);
    #elif defined(BOOST_HAS_GETTIMEOFDAY)
        struct timeval tv;
        gettimeofday(&tv, 0);
        xtp.sec = tv.tv_sec;
        xtp.nsec = tv.tv_usec * 1000;
    #elif defined(BOOST_HAS_CLOCK_GETTIME)
        timespec ts;
        clock_gettime(CLOCK_REALTIME, &ts);
        xtp.sec = ts.tv_sec;
        xtp.nsec = ts.tv_nsec;
    #else
    #error "time::get() implementation for your platform undefined"
    #endif
    }

    void    time::sleep(const time& xtp)
    {
        time    t,dt;
        for (int foo=0; foo < 5; ++foo)
        {
            time::get(t);
            if(!(t<xtp))
                return;
            dt  =   xtp-t;
        //FIXME: this is not quite clean. nicer way to detect if we are on windows
        //       but for the moment its sufficient
            #if defined(BOOST_HAS_FTIME)
            //we get :-( only millisecond resolution
            ::Sleep(dt.milli_seconds());
            #elif defined(BOOST_HAS_NANOSLEEP)
            timespec ts ={static_cast<int>(dt.sec), static_cast<int>(dt.nsec)};
            nanosleep(&ts,0);
            #else //use usleep
            usleep(xtp.nsec/NSEC_PER_USEC+xtp.sec*USEC_PER_SEC);
            #endif
        }
    }

    time::tick_type time::ticks()
    {
    static bool s_initialized=false;
    #if defined(BOOST_HAS_FTIME)
        /*FTIME version has too low resolution
          but maybe it is still useful

        FILETIME ft;
        static const boost::uint64_t TIMESPEC_TO_FILETIME_OFFSET =
            UINT64_C(116444736000000000);

        static boost::uint64_t s_init_time;

        if(!s_initialized)
        {
        #if defined(BOOST_NO_GETSYSTEMTIMEASFILETIME)
            {
                SYSTEMTIME st;
                GetSystemtime(&st);
                SystemtimeToFiletime(&st,&ft);
            }
        #else
            GetSystemtimeAsFiletime(&ft);
        #endif

            s_init_time =   ((static_cast<boost::uint64_t>(ft.dwHighDatetime) << 32)
                            + ft.dwLowDatetime)-TIMESPEC_TO_FILETIME_OFFSET;
            s_initialized=true;
        }

    #if defined(BOOST_NO_GETSYSTEMTIMEASFILETIME)
        {
            SYSTEMTIME st;
            GetSystemtime(&st);
            SystemtimeToFiletime(&st,&ft);
        }
    #else
        GetSystemtimeAsFiletime(&ft);
    #endif


        const boost::uint64_t ft64 = ((static_cast<boost::uint64_t>(ft.dwHighDatetime) << 32)
                                      + ft.dwLowDatetime)-TIMESPEC_TO_FILETIME_OFFSET;

        return (ft64 - s_init_time)/(NSEC_PER_MSEC/100);*/
        static LARGE_INTEGER s_init_time;
        static LARGE_INTEGER s_c_freq;
        if(!s_initialized)
        {
            QueryPerformanceCounter(&s_init_time);
            QueryPerformanceFrequency(&s_c_freq);
            //prepare this value so that division through it directly leads to milliseconds
            s_c_freq.QuadPart /= MSEC_PER_SEC;
            s_initialized = true;
        }
        LARGE_INTEGER t;
        QueryPerformanceCounter(&t);
        return static_cast<tick_type>((t.QuadPart - s_init_time.QuadPart)/s_c_freq.QuadPart);
    #elif defined(BOOST_HAS_GETTIMEOFDAY)
        static timeval s_init_time;
        if(!s_initialized)
        {
            #ifndef NDEBUG
            int res =
            #endif
            gettimeofday(&s_init_time, 0);
            assert(0 == res);
            assert(s_init_time.tv_sec >= 0);
            assert(s_init_time.tv_usec >= 0);
            s_initialized=true;
        }
        struct timeval tv;
        #ifndef NDEBUG
        int res =
        #endif
        gettimeofday(&tv, 0);
        assert(0 == res);
        assert(tv.tv_sec >= 0);
        assert(tv.tv_usec >= 0);
        return (tv.tv_sec-s_init_time.tv_sec)*MSEC_PER_SEC+ (tv.tv_usec-s_init_time.tv_usec)/USEC_PER_MSEC;
    #elif defined(BOOST_HAS_CLOCK_GETTIME)
        static timespec s_init_time;
        if(!s_initialized)
        {
            #ifndef NDEBUG
            int res =
            #endif
            clock_gettime(CLOCK_REALTIME,&s_init_time);
            assert(0 == res);
            s_initialized=true;
        }
        timespec ts;
        #ifndef NDEBUG
        int res =
        #endif
        clock_gettime(CLOCK_REALTIME, &ts);
        assert(0 == res);
        return (tv.tv_sec-s_init_time.tv_sec)*MSEC_PER_SEC+ (tv.tv_nsec-s_init_time.tv_nsec)/NSEC_PER_MSEC;
    #else
    #error "time::Ticks() implementation for your platform undefined"
    #endif
    }
} // namespace core
