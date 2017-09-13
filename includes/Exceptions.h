/* 
 * File:   Exceptions.h
 * Author: veraalva
 *
 * Created on March 24, 2016, 1:20 PM
 */

#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

namespace exceptions {

    /**
     * Exception to be thrown if there are problems opening files  
     */
    class FileHandledException : public std::exception {
    public:

        explicit FileHandledException(const char* message) : msg(message) {
        }

        explicit FileHandledException(const std::string& message) : msg(message) {
        }

        virtual ~FileHandledException() {
        }

        virtual const char* what() const throw () {
            return msg.c_str();
        }

    private:
        std::string msg;
    };

    /**
     * Exception to be thrown when an element is not found in the container
     */
    class NotFoundException : public std::exception {
    public:

        explicit NotFoundException(const char* message) : msg(message) {
        }

        explicit NotFoundException(const std::string& message) : msg(message) {
        }

        virtual ~NotFoundException() {
        }

        virtual const char* what() const throw () {
            return msg.c_str();
        }

    private:
        std::string msg;
    };
    
    class EmptyDatasetException : public std::exception {
    public:

        explicit EmptyDatasetException(const char* message) : msg(message) {
        }

        explicit EmptyDatasetException(const std::string& message) : msg(message) {
        }

        virtual ~EmptyDatasetException() {
        }

        virtual const char* what() const throw () {
            return msg.c_str();
        }

    private:
        std::string msg;
    };
    
    // ML_ERR_return_NAN
    class NANException : public std::exception {
    public:

        explicit NANException(const char* message) : msg(message) {
        }

        explicit NANException(const std::string& message) : msg(message) {
        }

        virtual ~NANException() {
        }

        virtual const char* what() const throw () {
            return msg.c_str();
        }

    private:
        std::string msg;
    };
}

#endif /* EXCEPTIONS_H */

