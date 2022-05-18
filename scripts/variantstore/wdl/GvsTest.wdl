import "./GvsUtil.wdl" as GvsUtil

workflow GvsTest {
    call GvsUtil.TerminateWorkflow {
        input: message = "Holy moly 'here' are a bunch of \"problematic\" quotes"
    }
}
