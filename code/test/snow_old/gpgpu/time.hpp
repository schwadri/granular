#ifndef _time_hpp_
#define _time_hpp_

/** \file time.hpp
 *  \author Adrian Schweizer
 *  \created  $Sa 18 Aug 06:32:47 pm CEST 2007 schwadri@SchwadriComp.local$
 *  \modified $So 19 Aug 05:21:30 pm CEST 2007 schwadri@SchwadriComp.local$
 *  \par FIXME
 *  - usage of integral types more consistent (all from boost or all from types.hpp)
 */

#include <iosfwd>
#include <boost/thread/xtime.hpp>
#include <boost/config.hpp>

namespace core {

    /** \brief Helper class to represent times and time intervals
     *
     *  This class supplies constructors and methods to represent an absolute point in time or a duration.
     *  The time is represented through two public member variables \a sec and \a nsec whose concrete type
     *  depend on the platform you are running on. Also depending on the platform absolute points in time
     *  may be defined starting from different starting points.
     *  All the comparision operators are implemented and times can also be added and substracted from each
     *  other. Be careful when doing substractions as the underlying data types to represent seconds and
     *  nanoseconds usually are unsigned.
     *
     *  \note time is intended to be used for high performance timing and not to represent times that lie back
     *  several years. Whereas on a posix platform this might still be possible on windows it is not because all
     *  times are relative to the system or process start time.
     *
     *  \par TODO
     *  - Usage Examples
     */
    class time //: public boost::xtime
    {
    public:
        typedef boost::xtime::xtime_sec_t       sec_type;   ///< the type used to represent seconds
        typedef boost::xtime::xtime_nsec_t      nsec_type;  ///< the type used to represend nanoseconds
        sec_type sec;
        nsec_type nsec;
        typedef int                             tick_type;  ///< the type used to represent "ticks" (milliseconds)
      
        static int const NSEC_PER_SEC;//   =   1000000000;
        static int const NSEC_PER_MSEC;//  =   1000000;
        static int const NSEC_PER_USEC;//  =   1000;
        static int const USEC_PER_MSEC;//  =   1000;
        static int const MSEC_PER_SEC;//   =   1000;
        static int const USEC_PER_SEC;//   =   1000000;

        enum sec_enum { Sec };
        enum milli_enum { Milli };
        enum micro_enum { Micro };
        enum nano_enum { Nano };

        /** \brief default c'tor
         *
         *  this does not initialize the members sec and nsec.
         */
        time()
        { }

        template<typename T>
            time(T seconds,sec_enum)
            {
                sec     =   static_cast<sec_type>(seconds);
                nsec    =   static_cast<nsec_type>( (seconds-sec)*NSEC_PER_SEC);
            }


        time(int seconds,sec_enum)
        {
            sec     =   static_cast<sec_type>(seconds);
            nsec    =   0;
        }

        template<typename T>
            time(T seconds,milli_enum)
            {
                sec     =   static_cast<sec_type>(seconds/MSEC_PER_SEC);
                nsec    =   static_cast<nsec_type>( (seconds-sec*MSEC_PER_SEC)*NSEC_PER_MSEC);
            }

        template<typename T>
            time(T seconds,micro_enum)
            {
                sec     =   static_cast<sec_type>(seconds/USEC_PER_SEC);
                nsec    =   static_cast<nsec_type>( (seconds-sec*USEC_PER_SEC)*NSEC_PER_USEC);
            }

        template<typename T>
            time(T seconds,nano_enum)
            {
                sec     =   static_cast<sec_type>(seconds/NSEC_PER_SEC);
                nsec    =   static_cast<nsec_type>(seconds-sec*NSEC_PER_SEC);
            }

        //millisecond is default
        template<typename T>
            explicit time(T milliseconds)
            {
                sec     =   static_cast<sec_type>(milliseconds / MSEC_PER_SEC);
                nsec    =   static_cast<nsec_type>((milliseconds - sec * MSEC_PER_SEC) * NSEC_PER_MSEC);
            }

        /*time(float seconds)
          {
          sec     =   static_cast<sec_type>(seconds);
          nsec    =   static_cast<nsec_type>( (seconds-sec)*NSEC_PER_SEC);
          }

          time(sec_type seconds)
          {
          sec     =   seconds;
          }*/

        time(sec_type seconds, nsec_type nano_seconds)
        {
            sec     =   seconds;// + nano_seconds/NSEC_PER_SEC ;
            nsec    =   nano_seconds;// % NSEC_PER_SEC;
        }

        time&   operator+=(const time& ot)
        {
            sec     +=  ot.sec;
            nsec    +=  ot.nsec;
            if(nsec >= NSEC_PER_SEC)
            {
                nsec    -=  NSEC_PER_SEC;
                ++sec;
            }
            else if(nsec <    0)
            {
                --sec;
                nsec    +=   NSEC_PER_SEC;
            }
            return *this;
        }

        time&   operator-=(const time& ot)
        {
            sec     -=  ot.sec;
            nsec    -=  ot.nsec;

            if(nsec <    0)
            {
                --sec;
                nsec    +=   NSEC_PER_SEC;
            }


            return *this;
        }

        tick_type milli_seconds() const
        {//FIXME:
			return static_cast<tick_type>(sec*MSEC_PER_SEC + nsec/NSEC_PER_MSEC);
        }

        int micro_seconds() const
        {//FIXME:
            return static_cast<tick_type>(sec*USEC_PER_SEC + nsec/NSEC_PER_USEC);
        }
        /** \brief return current time */
        static void get(time& xtp);
        /** \brief suspend the calling thread until \a xtp has been reached*/
        static void sleep(const time& xtp);
        static tick_type  ticks();    ///< returns current ticks (in milliseconds)

        #if defined(BOOST_HAS_CLOCK_GETTIME)
        inline void to_timespec(const time& xt, timespec& ts)
        {
            ts.tv_sec = static_cast<int>(xt.sec);
            ts.tv_nsec = static_cast<int>(xt.nsec);
            if(ts.tv_nsec >= time::NSEC_PER_SEC)
            {
                ts.tv_sec += ts.tv_nsec / time::NSEC_PER_SEC;
                ts.tv_nsec %= time::NSEC_PER_SEC;
            }
        }

        inline operator timespec()
        {
            timespec ts;
            ts.tv_sec = static_cast<int>(sec);
            ts.tv_nsec = static_cast<int>(nsec);
            if(ts.tv_nsec >= time::NSEC_PER_SEC)
            {
                ts.tv_sec += ts.tv_nsec / time::NSEC_PER_SEC;
                ts.tv_nsec %= time::NSEC_PER_SEC;
            }
            return ts;
        }
        #elif defined(BOOST_HAS_GETTIMEOFDAY)

        inline operator timeval()
        {
            timeval tv;
            tv.tv_sec = static_cast<int>(sec);
            tv.tv_usec = static_cast<int>(nsec/NSEC_PER_USEC);
            if(tv.tv_usec >=USEC_PER_SEC)
            {
                tv.tv_sec += tv.tv_usec / USEC_PER_SEC;
                tv.tv_usec %= USEC_PER_SEC;
            }
            return tv;
        }
        #endif

        static  const   time    PAST;           ///< a time Instance that represents a time far in the past
        static  const   time    FUTURE;         ///< a time Instance representing a time far into the future
    };

    time    operator+(const time& lt, const time& rt);

    time    operator-(const time& lt, const time& rt);

    bool    operator==(const time& lt, const time& rt);

    bool    operator!=(const time& lt, const time& rt);

    bool    operator<(const time& lt, const time& rt);

    bool    operator>(const time& lt, const time& rt);

    bool    operator<=(const time& lt, const time& rt);

    bool    operator>=(const time& lt, const time& rt);

    template<class _CharT, class _Traits>
        std::basic_ostream<_CharT, _Traits>&  operator <<(std::basic_ostream<_CharT, _Traits>& out,const time& t)
        {
            out << "[ sec:"<< t.sec << ", nsec:" << t.nsec << "]";
            return out;
        }

    /** \brief a little stopwatch which uses the time class
     *
     *  This class provides a minimal stopwatch to measure the elapsed time
     *  between calls to Start and Stop
     */
    class stop_watch
    {
        public:
            stop_watch()
            { }
            void start(const time& t)           { m_t_start = t;}
            const time&     start()
            {
                time::get(m_t_start);
                return m_t_start;
            }
            void            stop(const time& t) { m_t_stop = t; }
            const time&     stop()
            {
                time::get(m_t_stop);
                return m_t_stop;
            }
            time            duration()          { return m_t_stop - m_t_start; }
        private:
            time    m_t_start;
            time    m_t_stop;
    };

}   // namespace core


#endif // _time_hpp_

