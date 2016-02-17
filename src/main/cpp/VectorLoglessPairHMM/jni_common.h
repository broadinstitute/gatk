#ifndef JNI_COMMON_H
#define JNI_COMMON_H

#include <jni.h>

//#define ENABLE_ASSERTIONS 1
//#define DEBUG 1
//#define DUMP_TO_SANDBOX 1

#define GET_BYTE_ARRAY_ELEMENTS env->GetByteArrayElements
#define RELEASE_BYTE_ARRAY_ELEMENTS env->ReleaseByteArrayElements
#define JNI_RO_RELEASE_MODE JNI_ABORT
#define GET_DOUBLE_ARRAY_ELEMENTS env->GetDoubleArrayElements
#define RELEASE_DOUBLE_ARRAY_ELEMENTS env->ReleaseDoubleArrayElements

#endif  //ifndef JNI_COMMON_H
