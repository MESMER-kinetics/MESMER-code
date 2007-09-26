if HAVE_FORTRAN
AM_CPPFLAGS = -I$(top_builddir) -I$(top_builddir)/include -I$(top_srcdir)/include
LDADD=libqdmod.a $(top_builddir)/src/libqd.a $(FCLIBS)

if UPCASE_MODULE
QDMOD=QDMODULE.$(module_ext)
DDMOD=DDMODULE.$(module_ext)
else
QDMOD=qdmodule.$(module_ext)
DDMOD=ddmodule.$(module_ext)
endif

lib_LIBRARIES = libqdmod.a
libqdmod_a_SOURCES = ddmod.f qdmod.f f_dd.cpp f_qd.cpp main.cpp
qdmod.o: ddmod.o
$(QDMOD): qdmod.o
$(DDMOD): ddmod.o

pkglib_DATA = $(QDMOD) $(DDMOD)

DEMO=quaderq quadgsq2d quadgsq quadtsq quadtsq2d
TIMER = dd_timer qd_timer
EXTRA_PROGRAMS=$(DEMO) $(TIMER)

quaderq_SOURCES = tquaderq.f
nodist_quaderq_SOURCES=second.f
quaderq_LINK=$(CXXLINK)

quadgsq2d_SOURCES = tquadgsq2d.f
nodist_quadgsq2d_SOURCES=second.f
quadgsq2d_LINK=$(CXXLINK)

quadgsq_SOURCES = tquadgsq.f
nodist_quadgsq_SOURCES=second.f
quadgsq_LINK=$(CXXLINK)

quadtsq_SOURCES = tquadtsq.f
nodist_quadtsq_SOURCES=second.f
quadtsq_LINK=$(CXXLINK)

quadtsq2d_SOURCES = tquadtsq2d.f
nodist_quadtsq2d_SOURCES=second.f
quadtsq2d_LINK=$(CXXLINK)

dd_timer_SOURCES = dd_timer.f
nodist_dd_timer_SOURCES = second.f
dd_timer_LINK=$(CXXLINK)

qd_timer_SOURCES = qd_timer.f
nodist_qd_timer_SOURCES = second.f
qd_timer_LINK=$(CXXLINK)

demo: $(DEMO)

CLEANFILES = $(EXTRA_PROGRAMS) *.mod

endif

MAINTAINERCLEANFILES = dd_timer.f qd_timer.f
EXTRA_DIST = dd_timer.f qd_timer.f zz_timer.f Makefile.sample

dd_timer.f: zz_timer.f
	sed -e 's/@@/dd/' zz_timer.f >dd_timer.f

qd_timer.f: zz_timer.f
	sed -e 's/@@/qd/' zz_timer.f >qd_timer.f

if HAVE_FORTRAN
.PHONY: demo
endif

