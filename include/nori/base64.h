#if !defined(__NORI_BASE64_H)
#define __NORI_BASE64_H

#include <string>
#include <nori/common.h>

NORI_NAMESPACE_BEGIN

extern std::string base64_encode(unsigned char const* , unsigned int len);
extern std::string base64_decode(std::string const& s);

NORI_NAMESPACE_END

#endif

