/***************************************************************************
 * Copyright (C) 2023 Francesco Florian
 * This file is part of SkeletAlg.
 *
 * SkeletAlg is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SkeletAlg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SkeletAlg.  If not, see <http://www.gnu.org/licenses/>.
 *
 * The copyright holders give you permission to combine SkeletAlg
 * with code included in the standard release of Netgen (from Joachim
 * Schöberl), METIS (from George Karypis at the University of
 * Minnesota), OpenCASCADE (from Open CASCADE S.A.S) and ParaView
 * (from Kitware, Inc.) under their respective licenses. You may copy
 * and distribute such a system following the terms of the GNU GPL for
 * Gmsh and the licenses of the other code concerned, provided that
 * you include the source code of that other code when and as the GNU
 * GPL requires distribution of source code.
 *
 * Note that people who make modified versions of SkeletAlg are not
 * obligated to grant this special exception for their modified
 * versions; it is their choice whether to do so. The GNU General
 * Public License gives permission to release a modified version
 * without this exception; this exception also makes it possible to
 * release a modified version which carries forward this exception.
 *
 * Additional permission under GNU GPL version 3 section 7
 * If you modify this Program, or any covered work, by linking or combining it
 * with H2Lib (https://github.com/H2Lib/H2Lib), (or a modified version of that
 * library), containing parts covered by the terms all right reserved,the
 * licensors of this Program grant you additional permission to convey the
 * resulting work.
 ***************************************************************************/
#ifndef TIMING_H
#define TIMING_H

#include <cmath>
#include <ctime>
#include <list>
#include <string>

typedef struct _stopwatch stopwatch;

namespace MainNamespace {
    /**
     * Get elapsed time info (either with monotonic wall clock or cpu time)
     */
    class BaseTimer {
    public:
        virtual ~BaseTimer() {}         //!< \brief Destructor. One is needed to be inherited.
        virtual float read();           //!< \brief Read the elapsed time at the last stop().
        virtual void reset();           //!< \brief Reset the timer.
        virtual void resume() = 0;      //!< \brief Resume the timer.
        virtual void start();           //!< \brief Reset and start the timer.
        virtual void stop() = 0;        //!< \brief Store the current time in the timer.
    protected:
        std::clock_t total_ = 0;        //!< Sum all measured times since reset.
    };

    /**
     * Wall clock timer wrapper for the H2lib stopwatch
     */
    class BoermTimer : public BaseTimer {
    public:
        BoermTimer();                                   //! \brief New wall-clock timer
        BoermTimer(BoermTimer const &other) = default;  //!< \brief Explicitly enable default copy constructor.
        virtual ~BoermTimer() override; //! \brief destroy the timer data.
        virtual void resume() override;
        virtual void stop() override;   
    private:
        stopwatch *timer_;      //!< Pointer to the c struct holding the timer.
    };
  
    /**
     * Timer measuring cpu time
     */
    class CpuTimer : public BaseTimer {
    public:
        virtual ~CpuTimer() override = default;     //!< \brief Explicitly enable default destructor.
        virtual void resume() override;
        virtual void stop() override;
        virtual float read() override;
    private:
        std::clock_t start_ = 0;//!< Time when the timer was started.
    };

    /**
     * Collection of timers measuring the same interval in (possibly) different ways
     */
    class Clocks {
    private:
        std::list<BaseTimer*> timers_;
        std::string format_;
    public:
        ~Clocks();
        void format(std::string formatstring);
        void addTimer(BaseTimer* timer);
        void resume();
        void start();
        void stop();
        std::string readAll();
    };
}

#endif /* TIMING_H */
