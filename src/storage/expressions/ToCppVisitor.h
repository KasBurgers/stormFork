#pragma once

#include <sstream>
#include <unordered_map>

#include "src/storage/expressions/Variable.h"
#include "src/storage/expressions/ExpressionVisitor.h"

namespace storm {
    namespace expressions {
        class Expression;
        
        class ToCppTranslationOptions {
        public:
            ToCppTranslationOptions(std::string const& prefix = "", std::string const& valueTypeCast = "");
            
            std::string const& getPrefix() const;
            
            void setSpecificPrefixes(std::unordered_map<storm::expressions::Variable, std::string> const& prefixes);
            std::unordered_map<storm::expressions::Variable, std::string> const& getSpecificPrefixes() const;
            void clearSpecificPrefixes();
            
            bool hasValueTypeCast() const;
            std::string const& getValueTypeCast() const;
            void clearValueTypeCast();
            
        private:
            std::string valueTypeCast;
            std::string prefix;
            std::unordered_map<storm::expressions::Variable, std::string> specificPrefixes;
        };
        
        class ToCppVisitor : public ExpressionVisitor {
        public:
            std::string translate(storm::expressions::Expression const& expression, ToCppTranslationOptions const& options = ToCppTranslationOptions());
            
            virtual boost::any visit(IfThenElseExpression const& expression, boost::any const& data) override;
            virtual boost::any visit(BinaryBooleanFunctionExpression const& expression, boost::any const& data) override;
            virtual boost::any visit(BinaryNumericalFunctionExpression const& expression, boost::any const& data) override;
            virtual boost::any visit(BinaryRelationExpression const& expression, boost::any const& data) override;
            virtual boost::any visit(VariableExpression const& expression, boost::any const& data) override;
            virtual boost::any visit(UnaryBooleanFunctionExpression const& expression, boost::any const& data) override;
            virtual boost::any visit(UnaryNumericalFunctionExpression const& expression, boost::any const& data) override;
            virtual boost::any visit(BooleanLiteralExpression const& expression, boost::any const& data) override;
            virtual boost::any visit(IntegerLiteralExpression const& expression, boost::any const& data) override;
            virtual boost::any visit(RationalLiteralExpression const& expression, boost::any const& data) override;
            
        private:
            std::stringstream stream;
        };
        
    }
}
