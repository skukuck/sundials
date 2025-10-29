# TODO(CJB): Since litgen is GPLv3, this script might have to be GPLv3.
# Will need to determine if this is the case or not.
# The outputs of the script, i.e. the generated code, are definitely
# not subject to GPLv3 though, and can use our standard license.

import copy
import re
from typing import Optional
from codemanip import code_utils
from srcmlcpp.cpp_types import CppParameter
from litgen.internal.adapt_function_params._lambda_adapter import LambdaAdapter
from litgen.internal.adapted_types import AdaptedFunction, AdaptedParameter
from .utils import is_array_param


def adapt_sundials_types_returns_to_shared_ptr(
    adapted_function: AdaptedFunction,
) -> Optional[LambdaAdapter]:
    """
    Adapts functions that take pointer_type* (pointer to pointer_type) parameters to instead accept
    std::vector<pointer_type> in the Python interface. Handles multiple such parameters in one function.
    """
    options = adapted_function.options
    sundials_pointer_types = options.sundials_pointer_types
    function_name = adapted_function.cpp_adapted_function.function_name
    return_type = adapted_function.cpp_adapted_function.return_type

    print(f"processing {function_name}")
    print(f"    return type: {return_type}")

    needs_adapt = False
    is_tuple = False

    if return_type.str_return_type().startswith("std::tuple"):
        # Check if any std::tuple return type argument needs adaptation
        tuple_args_match = re.match(r"std::tuple\s*<(.+?)>", return_type.str_return_type())
        needs_adapt = (
            any(
                arg.strip() in sundials_pointer_types
                for arg in tuple_args_match.group(1).split(",")
            )
            if tuple_args_match
            else False
        )
        is_tuple = True
    elif return_type.str_return_type() in sundials_pointer_types:
        needs_adapt = True

    if not needs_adapt:
        return None

    # Ensure SUNContext is kept alive while this object is alive
    # The easiest way to do this without modifying litgen is to inject the
    # keep_alive statement as a comment.
    for idx, param in enumerate(adapted_function.adapted_parameters()):
        if "SUNContext" in param.cpp_element().decl.cpp_type.typenames:
            adapted_function.cpp_element().cpp_element_comments.add_eol_comment(
                f"nb::keep_alive<0, {idx+1}>()"
            )

    # Wrap return type in shared_ptr
    lambda_adapter = LambdaAdapter()
    lambda_adapter.new_function_infos = copy.deepcopy(adapted_function.cpp_adapted_function)
    old_function_params: list[AdaptedParameter] = adapted_function.adapted_parameters()
    new_function_params: list[CppParameter] = old_function_params
    old_return_type = adapted_function.cpp_adapted_function.return_type

    if is_tuple:
        # new_return_type = copy.deepcopy(old_return_type)
        # new_return_type.modifiers = []
        # new_return_type.specifiers = []
        # new_return_type.typenames = [f"std::shared_ptr<std::remove_pointer_t<{old_return_type}>>"]
        # lambda_adapter.new_function_infos.return_type = new_return_type

        # # Note: litgen always generates the statement `return lambda_result;` if the function return type is not void
        # # and places it after lambda_output_code. As such, this is going to result in two return statements. First,
        # # the one beloe that we want and then the one litgen always adds. This might cause compiler warnings, but is harmless
        # # otherwise.
        # lambda_output_code = f"""
        #     return our_make_shared<std::remove_pointer_t<{str(old_return_type)}>, {str(old_return_type)}Deleter>(lambda_result);
        # """

        # lambda_adapter.lambda_output_code += (
        #     code_utils.unindent_code(lambda_output_code, flag_strip_empty_lines=True) + "\n"
        # )

        # for param in new_function_params:
        #     lambda_adapter.adapted_cpp_parameter_list.append(
        #         f"{param.cpp_element().decl.decl_name}"
        #     )

        # lambda_adapter.lambda_name = (
        #     adapted_function.cpp_adapted_function.function_name
        #     + f"_adapt_return_type_to_shared_ptr"
        # )
        return None
    else:
        new_return_type = copy.deepcopy(old_return_type)
        new_return_type.modifiers = []
        new_return_type.specifiers = []
        new_return_type.typenames = [f"std::shared_ptr<std::remove_pointer_t<{old_return_type}>>"]
        lambda_adapter.new_function_infos.return_type = new_return_type

        # Note: litgen always generates the statement `return lambda_result;` if the function return type is not void
        # and places it after lambda_output_code. As such, this is going to result in two return statements. First,
        # the one beloe that we want and then the one litgen always adds. This might cause compiler warnings, but is harmless
        # otherwise.
        lambda_output_code = f"""
            return our_make_shared<std::remove_pointer_t<{str(old_return_type)}>, {str(old_return_type)}Deleter>(lambda_result);
        """

        lambda_adapter.lambda_output_code += (
            code_utils.unindent_code(lambda_output_code, flag_strip_empty_lines=True) + "\n"
        )

        for param in new_function_params:
            lambda_adapter.adapted_cpp_parameter_list.append(
                f"{param.cpp_element().decl.decl_name}"
            )

        lambda_adapter.lambda_name = (
            adapted_function.cpp_adapted_function.function_name
            + f"_adapt_return_type_to_shared_ptr"
        )

    return lambda_adapter
