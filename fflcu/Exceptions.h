/*
 * Exceptions.h
 *
 *  Created on: Jan 31, 2012
 *      Author: alexander
 */

#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_
#include <string>
	struct WrongDeviceNumberException {
		WrongDeviceNumberException (int n) {
			deviceNumber = n;
		}
		unsigned int deviceNumber;
	};

	struct TooDenseSystemException{
		TooDenseSystemException(unsigned int lastPart){
			lastParticleNumber = lastPart;
		}
		unsigned int lastParticleNumber;
	};

	struct DeviceMemoryException{
		DeviceMemoryException(std::string message){
			this->message = message;
		}
		std::string message;
	};

	struct HostMemoryAllocationException{
		HostMemoryAllocationException(std::string message){
			this->message = message;
		}
		std::string message;
	};

	struct DeviceMemoryAllocationException: public DeviceMemoryException{
		DeviceMemoryAllocationException(std::string message):DeviceMemoryException(message){}
	};

	struct DeviceMemoryCopyException: public DeviceMemoryException{
		DeviceMemoryCopyException(std::string message):DeviceMemoryException(message){}
	};

	struct DeviceMemCpyToSymbolException: public DeviceMemoryException{
		DeviceMemCpyToSymbolException(std::string message):DeviceMemoryException(message){}
	};

	struct DeviceMemsetException: public DeviceMemoryException{
		DeviceMemsetException(std::string message):DeviceMemoryException(message){}
	};

#endif /* EXCEPTIONS_H_ */
