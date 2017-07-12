 /**
  * The MIT License (MIT)
  *
  * Copyright (c) 2014 Behrouz Babaki, Tias Guns, Siegfried Nijssen
  *
  * Permission is hereby granted, free of charge, to any person obtaining a copy
  * of this software and associated documentation files (the "Software"), to deal
  * in the Software without restriction, including without limitation the rights
  * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  * copies of the Software, and to permit persons to whom the Software is
  * furnished to do so, subject to the following conditions:
  *
  * The above copyright notice and this permission notice shall be included in
  * all copies or substantial portions of the Software.
  *
  * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  * THE SOFTWARE.
  */

#ifndef WRAPPERS_H
#define WRAPPERS_H

#include <scip/def.h>

#if SCIP_VERSION == 301
#define WrapperSCIPgetLPBranchCands(scip,lpcands,lpcandssol,lpcandsfrac,nlpcands,npriolpcands) SCIPgetLPBranchCands(scip,lpcands,lpcandssol,lpcandsfrac,nlpcands,npriolpcands) 
#elif SCIP_VERSION >= 310
#define WrapperSCIPgetLPBranchCands(scip,lpcands,lpcandssol,lpcandsfrac,nlpcands,npriolpcands) SCIPgetLPBranchCands(scip,lpcands,lpcandssol,lpcandsfrac,nlpcands,npriolpcands,NULL) 
#endif

#if SCIP_VERSION < 320
#define WrapperObjConshdlr(scip, name, desc, sepapriority, enfopriority, checkpriority, sepafreq, propfreq, eagerfreq, maxprerounds, delaysepa, delayprop, delaypresol, needscons, proptiming, presoltiming) ObjConshdlr(scip, name, desc, sepapriority, enfopriority, checkpriority, sepafreq, propfreq, eagerfreq, maxprerounds, delaysepa, delayprop, delaypresol, needscons, proptiming)
#else
#define WrapperObjConshdlr(scip, name, desc, sepapriority, enfopriority, checkpriority, sepafreq, propfreq, eagerfreq, maxprerounds, delaysepa, delayprop, delaypresol, needscons, proptiming, presoltiming) ObjConshdlr(scip, name, desc, sepapriority, enfopriority, checkpriority, sepafreq, propfreq, eagerfreq, maxprerounds, delaysepa, delayprop, needscons, proptiming, presoltiming)
#endif

#endif
